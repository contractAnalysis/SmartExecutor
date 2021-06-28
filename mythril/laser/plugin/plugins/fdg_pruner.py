from copy import copy

from virtualenv.config.convert import NoneType


from fdg.FDG_3d_array import FDG
from fdg.funtion_info import Function_info
from fdg.sequence import Sequence
from mythril.laser.ethereum.svm import LaserEVM
from mythril.laser.plugin.interface import LaserPlugin
from mythril.laser.plugin.builder import PluginBuilder
from mythril.laser.plugin.plugins.dependency_pruner import get_dependency_annotation
from mythril.laser.plugin.signals import PluginSkipState

from mythril.laser.ethereum.state.global_state import GlobalState
from mythril.laser.ethereum.transaction.transaction_models import (
    ContractCreationTransaction,
)


import logging
import fdg.FDG_global
import time
import numpy as np
log = logging.getLogger(__name__)


class FDG_prunerBuilder(PluginBuilder):
    name = "fdg-pruner"

    def __call__(self, *args, **kwargs):
        return FDG_pruner()


class FDG_pruner(LaserPlugin):
    """ """

    def __init__(self):
        """Creates FDG pruner"""
        self._reset()


    def _reset(self):
        self._iteration_ = 0
        self._depth_=0

        self._level_=0
        self._level_control_list=[]
        self._level_control_list.append(0)
        self._level_control_index=0

        self.solidity=''
        self.contract=''
        self.FDG=None
        self.seq_object=None
        self.ftn_start_nodes=[] # save the functions reached at the depth 1
        self.ftn_wr_ftns_pc={} # save the pc of instruction for matching each function

        self.flag_changing_pc=False
        self.ftn_pc={}
        self.ftn_target={}
        self.ftn_covered_mark=None
        self.OS_states={}
        self.ftn_pairs={}
        self.info_1st={}
        self.info_following={}

    def initialize(self, symbolic_vm: LaserEVM) -> None:
        """Initializes the FDG_pruner
        :param symbolic_vm
        """
        self._reset()

        @symbolic_vm.laser_hook("start_sym_exec")
        def start_sym_exec_hook():

            # initialize FDG
            self.solidity=fdg.FDG_global.solidity_path
            self.contract=fdg.FDG_global.contract
            function_info = Function_info(self.solidity, self.contract)
            self.FDG= FDG(function_info.functions_dict_slither())

            self.ftn_covered_mark=[0]*self.FDG.num_ftn # to record which functions are covered
            self.ftn_covered_mark[0]=1 # constructor is covered

        @symbolic_vm.laser_hook("stop_sym_exec")
        def stop_sym_exec_hook():
            print(f' End of symbolic execution')
            ftn_not_covered=list(np.where(np.array(self.ftn_covered_mark) == 0)[0])
            if ftn_not_covered:
                ftns=[self.FDG.index_to_ftn[ftn_idx] for ftn_idx in ftn_not_covered]
                print(f'Functions not covered={ftns}')



        #-------------------------------------------------
        ''' 
          new hook methods for changing laserEVM instance
        - add states 
        - save states        
        '''
        #-------------------------------------------------
        @symbolic_vm.laser_hook("start_sym_trans_laserEVM")
        def start_sym_trans_hook_laserEVM(laserEVM:LaserEVM):
            """
            add states to laserEVM.open_states so that they can be used
            as base states in the next iteration of symbolic transaction
            :param laserEVM: instance of LaserEVM
            :return:
            """

            # generate Sequence object
            if self._iteration_ == self.FDG.depth_all_ftns_reached + 1:
                # sequence genration
                # build new fdg
                self.FDG.build_fdg_3d_array_new(self.ftn_pairs)

                ftn_idx_not_covered = list(np.where(np.array(self.ftn_covered_mark) == 0)[0])
                # create an Sequence object
                seq_object = Sequence(self.FDG, ftn_idx_not_covered, 5)

                # get ready for sequence generation(parpare some data for sequence generation)
                seq_object.prepare_sequence_generation()

                self.seq_object=seq_object

                self._level_control_list.append(self._level_+1)

            if self._iteration_>self.FDG.depth_all_ftns_reached:
                self._level_control_index+=1

                if len(self._level_control_list)==self._level_control_index:
                    print(f'end of all levels')
                    fdg.FDG_global.transaction_count=self._iteration_
                    return
                if self._level_control_list[self._level_control_index]>0:
                    # the start of a new level
                    # get the level
                    self._level_=self._level_control_list[self._level_control_index]
                    # update functions not covered
                    self.seq_object.ftn_idx_not_covered=list(np.where(np.array(self.ftn_covered_mark) == 0)[0])
                    # get sequences for functions not covered
                    sequences_dict = self.seq_object.get_sequences_by_level(self._level_)
                    print(f'sequences_dict ={sequences_dict }')

                    info_1st_temp,info_following_temp=self.seq_object.organize_sequences_for_execution_by_level(sequences_dict,self._level_)
                    if info_1st_temp:
                        states_for_seq_execution =[]
                        for depth,ftn_dict in info_1st_temp.items():
                            for key_ftn in ftn_dict.keys():
                                # the states at depth_all_ftns_reached are already there.
                                if depth-1 in self.OS_states.keys():
                                    if key_ftn in self.OS_states[depth-1].keys():
                                        states_for_seq_execution += self.OS_states[depth - 1][key_ftn]


                        # place states to laserEVM instance
                        if len(states_for_seq_execution)>0:
                            laserEVM.open_states = states_for_seq_execution

                            # replace functions to be executd with its associated PCs
                            for key,values in info_1st_temp.items():
                                for k,ftn_list in values.items():
                                    pc_values=[]
                                    for ftn_idx in ftn_list:
                                        ftn_name = self.FDG.index_to_ftn[ftn_idx]
                                        ftn_4byte = self.FDG.ftn_to_selector[ftn_name]
                                        pc_values.append(self.ftn_pc[ftn_4byte])
                                    pc_values.sort()
                                    info_1st_temp[key][k]=pc_values
                            self.info_1st[self._level_]=info_1st_temp

                            if info_following_temp:
                                for key,values in info_following_temp.items():
                                    self._level_control_list.append(-key) #  indicate the case to continue sequence execution
                                    for k,ftn_list in values.items():
                                        pc_values=[]
                                        for ftn_idx in ftn_list:
                                            ftn_name = self.FDG.index_to_ftn[ftn_idx]
                                            ftn_4byte = self.FDG.ftn_to_selector[ftn_name]
                                            pc_values.append(self.ftn_pc[ftn_4byte])
                                        pc_values.sort()
                                        info_following_temp[key][k] = pc_values

                                self.info_following[self._level_]=info_following_temp

                            else: self._level_control_list.append(self._level_ + 1)

                        else:
                            # go the next level
                            self._level_control_list.append(self._level_ + 1)
                            laserEVM.open_states = []
                            print(f'no states for sequence execution at level: {self._level_}')
                    else:
                        # go the next level
                        self._level_control_list.append(self._level_+1)
                        laserEVM.open_states=[]
                else:
                    # do not do anything here, just continue sequence execution
                    pass


        @symbolic_vm.laser_hook("stop_sym_trans_laserEVM")
        def stop_sym_trans_hook_laserEVM(laserEVM:LaserEVM):
            """
            - save states at the end of symbolic transactions for sequence execution
            - only need to states from depth 1 to self.FDG.depth_all_ftns_reached -1
            - the states at self.FDG.depth_all_ftns_reached are already saved in laserEVM
            - some saved states are used as initial states in sequence execution

            :param laserEVM:
            :return:
            """
            # assume that states are annotated by function sequences from which they are generated

            # get functions reached from constructor at depth 1
            if self._depth_ == 1:
                # only consider functions in the FDG
                for state in laserEVM.open_states:
                    if not state.constraints.is_possible: continue
                    if state.node.function_name in self.FDG.ftn_to_index.keys():
                        ftn_idx = self.FDG.ftn_to_index[state.node.function_name]
                        self.ftn_start_nodes.append(ftn_idx)
                        self.ftn_covered_mark[ftn_idx] = 1
                        # save states at depth 1
                        if ftn_idx not in self.OS_states[self._depth_ - 1].keys():
                            self.OS_states[self._depth_ - 1][ftn_idx] = [copy(state)]
                        else:
                            self.OS_states[self._depth_- 1][ftn_idx] += [copy(state)]


            # save states from depth 2 to self._depth__all_ftns_reached
            if self._depth_<=self.FDG.depth_all_ftns_reached and self._depth_>1:
                    for state in laserEVM.open_states:
                        if not state.constraints.is_possible:continue
                        if state.node.function_name in self.FDG.ftn_to_index.keys():
                            ftn_name=state.node.function_name
                            ftn_idx=self.FDG.ftn_to_index[ftn_name]

                            if ftn_idx not in self.OS_states[self._depth_ - 1].keys():
                                self.OS_states[self._depth_ - 1][ftn_idx] = [copy(state)]
                            else:
                                self.OS_states[self._depth_ - 1][ftn_idx] += [copy(state)]

                            self.ftn_covered_mark[ftn_idx] = 1  # record the function covered


            # check if target functions are covered or not
            if self._iteration_ > self.FDG.depth_all_ftns_reached:
                for state in laserEVM.open_states:
                    if not state.constraints.is_possible: continue
                    if state.node.function_name in self.FDG.ftn_to_index.keys():
                        ftn_idx = self.FDG.ftn_to_index[state.node.function_name]
                        self.ftn_covered_mark[ftn_idx] = 1

        @symbolic_vm.laser_hook("start_sym_trans")
        def start_sym_trans_hook():

            self._iteration_+=1
            self._depth_ += 1

            # define variables to save states and function pairs
            if self._depth_<=self.FDG.depth_all_ftns_reached and self._depth_>=1:
                self.OS_states[self._depth_-1]={}
                if self._depth_>=2:
                    self.ftn_pairs[self._depth_-1]=[]

        @symbolic_vm.laser_hook('stop_sym_trans')
        def stop_sym_trans_hook():
            '''
            build FDG after the first symbolic transaction
            because which functions are reachable from constructor are needed to build FDG
            as a way to prune false dependency.
            '''
            # at depth 1, build FDG
            if self._depth_==1:
                self.ftn_start_nodes=list(set(self.ftn_start_nodes))
                if len(self.ftn_start_nodes)==0:
                    fdg.FDG_global.transaction_count = self._iteration_
                    print(f'No function reachable from constructor')
                    return

                # build FDG
                self.FDG.build_fdg_3d_array(self.ftn_start_nodes)

                # update the value fdg.FDG_global.depth_all_ftns_reached
                fdg.FDG_global.depth_all_ftns_reached=self.FDG.depth_all_ftns_reached


                # replace dependent function indices with their PCs,
                # which are entry points in 'function dispatch' to evaluate functions.
                for key, values in self.FDG.graph_forward.items():
                    if isinstance(values, list):
                        if len(values) > 0:
                            pc_values = []
                            for ftn_idx in values:
                                ftn_name = self.FDG.index_to_ftn[ftn_idx]
                                ftn_4byte = self.FDG.ftn_to_selector[ftn_name]
                                if ftn_4byte in self.ftn_pc:
                                    pc_values.append(self.ftn_pc[ftn_4byte])
                                else:
                                    print(f'function {ftn_name} can not be found in {self.contract}')
                            pc_values.sort()
                            self.ftn_wr_ftns_pc[key] = pc_values


            # control when to stop by modifying transaction_count
            if self.FDG.depth_all_ftns_reached==1:
                fdg.FDG_global.transaction_count = self._iteration_
                print(f'FDG has depth of 1')
            if self._iteration_>=self.FDG.depth_all_ftns_reached:
                # stop symbolic execution if all functions are covered
                if self.ftn_covered_mark.count(0)==0:
                    # set to the current iteration, so that execution engine can stop
                    fdg.FDG_global.transaction_count = self._iteration_




        @symbolic_vm.laser_hook("execute_state")
        def execute_state_hook(global_state: GlobalState):
            '''
            do some thing on each state
            :param global_state:
            :return:
            '''
            pass

        #-------------------------------------------------
        ''' 
        - handle everything related to driving symbolic execution engine to 
        points required to evaluate specified functions
        - basically, this is done by:
            saving PCs for evaluating functions
            and changing machine state PC to PCs associated to specified functions
        '''
        #-------------------------------------------------
        @symbolic_vm.pre_hook("PUSH4")
        def push4_hook(state: GlobalState):
            # get the pc of opcode starting matching for each function
            if state.environment.active_function_name == 'fallback' \
            and self._depth_ == 1:
                ftn_selector = state.instruction['argument'][2:]
                # assert (ftn_selector not in self.ftn_pc.keys())
                if ftn_selector not in self.ftn_pc.keys():
                    self.ftn_pc[ftn_selector]=state.mstate.pc



        @symbolic_vm.post_hook("DUP1")
        def dup1_hook(state: GlobalState):
            if state.environment.active_function_name == 'fallback'\
            and self._iteration_ >= 2:
                # only consider DUP1 within a specified range
                pc_here = state.mstate.pc
                if pc_here < self.ftn_target['pc_interval_start']: return
                if pc_here > self.ftn_target['mismatch']: return


                annotations = get_dependency_annotation(state)
                ftn_seq=annotations.ftn_seq
                pre_ftn_name=ftn_seq[-1] # get the function name
                if pre_ftn_name not in self.FDG.ftn_to_index.keys():
                    state.mstate.pc=self.ftn_target['mismatch']
                    return
                pre_ftn_name_idx=self.FDG.ftn_to_index[pre_ftn_name]

                if self._iteration_<=self.FDG.depth_all_ftns_reached:
                    # for  FDG guided execution process
                    if pre_ftn_name_idx in self.ftn_wr_ftns_pc.keys():
                        state.mstate.pc=self.ftn_wr_ftns_pc[pre_ftn_name_idx][0]
                    else:
                        state.mstate.pc = self.ftn_target['mismatch']

                else:
                    # for sequence execution
                    # for the first step
                    if self._level_==self._level_control_list[self._level_control_index]:
                        print(f'_level={self._level_}')
                        print(f'_state.mstate.pc={state.mstate.pc}')
                        if len(ftn_seq) in self.info_1st[self._level_].keys():
                            if pre_ftn_name_idx in self.info_1st[self._level_][len(ftn_seq)].keys():
                                state.mstate.pc = self.info_1st[self._level_][len(ftn_seq)][pre_ftn_name_idx][0]
                            else:pass
                        else:
                            state.mstate.pc = self.ftn_target['mismatch']
                        print(f'_state.mstate.pc={state.mstate.pc}')
                    else:
                        # get the step for sequence execution
                        seq_execution_step=-self._level_control_list[self._level_control_index]
                        if seq_execution_step in self.info_following[self._level_].keys():
                            if pre_ftn_name_idx in self.info_following[self._level_][seq_execution_step].keys():
                                state.mstate.pc = self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx][0]
                        else:
                            state.mstate.pc = self.ftn_target['mismatch']


        @symbolic_vm.post_hook("JUMPI")
        def jumpi_hook(state: GlobalState):
            if self._depth_>=2 and state.environment.active_function_name=='fallback':
                pc_here = state.mstate.pc
                if pc_here < self.ftn_target['pc_interval_start']: return
                if pc_here > self.ftn_target['mismatch']: return





                annotations = get_dependency_annotation(state)
                ftn_seq=annotations.ftn_seq

                pre_ftn_name = ftn_seq[-1]  # get the function name
                # if pre_ftn_name does not in FDG, then we are not interested in it, jump to block leading to revert.
                if pre_ftn_name not in self.FDG.ftn_to_index.keys():
                    state.mstate.pc=self.ftn_target['mismatch']
                    return
                pre_ftn_name_idx =self.FDG.ftn_to_index[pre_ftn_name]

                if self._iteration_<=self.FDG.depth_all_ftns_reached:
                    # for FDG guided execution process
                    if pre_ftn_name_idx in self.ftn_wr_ftns_pc.keys()\
                    and len(self.ftn_wr_ftns_pc[pre_ftn_name_idx])>=2:
                        for i in range(len(self.ftn_wr_ftns_pc[pre_ftn_name_idx])-1):
                            if state.mstate.pc < self.ftn_wr_ftns_pc[pre_ftn_name_idx][i + 1] and \
                                state.mstate.pc > self.ftn_wr_ftns_pc[pre_ftn_name_idx][i]:
                                state.mstate.pc = self.ftn_wr_ftns_pc[pre_ftn_name_idx][i + 1]
                                return
                        if state.mstate.pc < self.ftn_target['mismatch'] \
                            and state.mstate.pc > self.ftn_wr_ftns_pc[pre_ftn_name_idx][-1]:
                            state.mstate.pc = self.ftn_target['mismatch']
                    else:
                        if state.mstate.pc < self.ftn_target['mismatch']:
                            state.mstate.pc = self.ftn_target['mismatch']
                else:
                    # for sequence execution
                    if self._level_==self._level_control_list[self._level_control_index]:
                        if len(ftn_seq) in self.info_1st[self._level_].keys():

                            if pre_ftn_name_idx in self.info_1st[self._level_][len(ftn_seq)].keys():
                                if len(self.info_1st[self._level_][len(ftn_seq)][pre_ftn_name_idx])>=2:
                                    for i in range(len(self.info_1st[self._level_][len(ftn_seq)][pre_ftn_name_idx]) - 1):
                                        if state.mstate.pc < self.info_1st[self._level_][len(ftn_seq)][pre_ftn_name_idx][i + 1] and \
                                            state.mstate.pc > self.info_1st[self._level_][len(ftn_seq)][pre_ftn_name_idx][i]:
                                            state.mstate.pc = self.info_1st[self._level_][len(ftn_seq)][pre_ftn_name_idx][i + 1]
                                            return
                                    if state.mstate.pc < self.ftn_target['mismatch'] \
                                        and state.mstate.pc > self.info_1st[self._level_][len(ftn_seq)][pre_ftn_name_idx][-1]:
                                        state.mstate.pc = self.ftn_target['mismatch']

                                else:
                                    if state.mstate.pc < self.ftn_target['mismatch']:
                                        state.mstate.pc = self.ftn_target['mismatch']

                    else: # for sequence execution, not the fist step
                        # get the step for sequence execution
                        seq_execution_step=-self._level_control_list[self._level_control_index]
                        if seq_execution_step in self.info_following[self._level_].keys():
                            if pre_ftn_name_idx in self.info_following[self._level_][seq_execution_step].keys():
                                if len(self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx]) >= 2:
                                    for i in range(len(self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx]) - 1):
                                        if state.mstate.pc < self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx][
                                            i + 1] and \
                                            state.mstate.pc > self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx][i]:
                                            state.mstate.pc = self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx][i + 1]
                                            return
                                    if state.mstate.pc < self.ftn_target['mismatch'] \
                                    and state.mstate.pc > self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx][-1]:
                                        state.mstate.pc = self.ftn_target['mismatch']

                                else:
                                    if state.mstate.pc < self.ftn_target['mismatch']:
                                        state.mstate.pc = self.ftn_target['mismatch']

        @symbolic_vm.pre_hook("JUMPDEST")
        def jumpdest_hook(state: GlobalState):
            '''
            save the PC used when the function selector is not matched
            this PC leads to the block with revert as the final opcode
            :param state:
            :return:
            '''
            # get the pc of jumpdest, the start opcode of a block meaning the end of function mapping
            if state.environment.active_function_name == 'fallback':
                if self._depth_==1:
                    self.ftn_target['mismatch']=state.mstate.pc

        @symbolic_vm.pre_hook("CALLDATASIZE")
        def calldatasize_hook(state: GlobalState):
            '''
            save the PC used when the function selector is not matched
            this PC leads to the block with revert as the final opcode
            :param state:
            :return:
            '''
            # get the pc of jumpdest, the start opcode of a block meaning the end of function mapping
            if state.environment.active_function_name == 'fallback':
                if self._depth_==1:
                    self.ftn_target['pc_interval_start']=state.mstate.pc


        @symbolic_vm.pre_hook("STOP")
        def stop_hook(state: GlobalState):
            _transaction_end(state)

        @symbolic_vm.pre_hook("RETURN")
        def return_hook(state: GlobalState):
            _transaction_end(state)

        def _transaction_end(state: GlobalState) -> None:
            """
            - collect function pairs that the second function is executed with return or stop
            - function pair: (fa,fb), fb depends on fa. fa provides states as initial states for executing fb.
                if fb is executed without revert, fb depends on fa; otherwise, fb does not.

            :param state:
            """
            # get function pairs from depth 1 to self.FDG.depth_all_ftns_reached
            if self._depth_>=2 and self._depth_<=self.FDG.depth_all_ftns_reached:
                ftn_seq=get_dependency_annotation(state).ftn_seq
                ftn_idx_seq=[self.FDG.ftn_to_index[ftn_name] for ftn_name in ftn_seq[-2:] if ftn_name in self.FDG.ftn_to_index.keys() ]

                if len(ftn_idx_seq)<2:
                    # no need to add function pairs
                    # ( quite possible one of function is function not in FDG like fallback function)
                    return
                self.ftn_pairs[self._depth_-1]+=[(ftn_idx_seq[0],ftn_idx_seq[1])]



        @symbolic_vm.laser_hook("add_world_state")
        def world_state_filter_hook(state: GlobalState):
            if isinstance(state.current_transaction, ContractCreationTransaction):
                # Reset iteration variable
                self._iteration_ = 0
                self._depth_=0
                return



from copy import copy

from virtualenv.config.convert import NoneType

from fdg import utils
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
        self._depth_ = 0
        self._depth_=0

        self._level_=0
        self._level_control_list=[0] 
        self._level_control_index=0

        self.solidity=''
        self.contract=''
        self.FDG=None
        self.seq_object=None
        self.ftn_start_nodes=[] # save the functions reached at the depth 1
        self.fdg_pc={} # save the pc of instruction for matching each function

        self.selector_pc={}
        self.ftn_pc={}
        self.pc_ftn={}
        self.pc_control_interval={}
        self.gt_pc = []
        self.valid_pc_interval = []

        self.ftn_covered_mark=None
        self.ftn_not_covered_pc_dict={}
        self.ftn_wo_edges_not_covered=[]
        self.OS_states={}
        self.ftn_pairs={}

        self.info_1st={}
        self.info_following={}
        self.no_public_ftn=False
        self.all_ftn_pc_list=[]


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
            print(f'End of symbolic execution')         
            ftns=[self.FDG.index_to_ftn[ftn_i] for ftn_i, item in enumerate(self.ftn_covered_mark) if item==0]
            print(f'Functions not covered={ftns}')

        @symbolic_vm.laser_hook("start_sym_trans")
        def start_sym_trans_hook():
            self._depth_+=1

            # define variables to save states and function pairs
            if self._depth_ <= fdg.FDG_global.depth_all_ftns_reached+1:
                self.OS_states[self._depth_]={}
                if self._depth_>=2:
                    self.ftn_pairs[self._depth_-1]=[]

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


            # update and/or add data in self.fdg_pc for the next iteration in FDG exploration
            if self._depth_ <= fdg.FDG_global.depth_all_ftns_reached+1 and self._depth_>1:
                prt_no_child = []  # record funtions who do not have dependent functions
                if self._depth_-1==1:
                    prt_ftns=self.ftn_start_nodes
                else:
                    prt_ftns=self.OS_states[self._depth_-1].keys()


                for ftn_idx in prt_ftns:
                    child_nodes=[]
                    if ftn_idx in self.FDG.graph.keys():
                        child_nodes=self.FDG.graph[ftn_idx]
                    if len(child_nodes)==0:
                        prt_no_child.append(ftn_idx)
                        if self._depth_ == 2:  # only do this at depth 2
                            # if len(self.ftn_wo_edges_not_covered) > 0:
                            #     child_nodes_pc = [self.ftn_pc[ftn_i] for ftn_i in self.ftn_wo_edges_not_covered if
                            #                       ftn_i in self.ftn_pc.keys()]
                            #     self.fdg_pc[ftn_idx] = sorted(child_nodes_pc)
                            # else:
                            #     pc_covered = [pc for ftn_i, pc in self.ftn_pc.items() if
                            #                   self.ftn_covered_mark[ftn_i] == 1]
                            #     pc_list = [pc for pc in self.all_ftn_pc_list if pc not in pc_covered]
                            #     self.fdg_pc[ftn_idx] = sorted(pc_list)
                            pc_covered = [pc for ftn_i, pc in self.ftn_pc.items() if
                                          self.ftn_covered_mark[ftn_i] == 1]
                            pc_list = [pc for pc in self.all_ftn_pc_list if pc not in pc_covered]
                            if len(self.valid_pc_interval)>0:
                                pc_list=[pc for pc in pc_list if utils.pc_is_valid(pc,self.valid_pc_interval)]
                            self.fdg_pc[ftn_idx] = sorted(pc_list)

                        continue
                    if self._depth_<=fdg.FDG_global.depth_all_ftns_reached:
                        # update or add data in self.fdg_pc
                        if ftn_idx not in self.fdg_pc.keys() or len(self.ftn_wo_edges_not_covered)>0:
                            child_nodes += self.ftn_wo_edges_not_covered
                            child_nodes = list(set(child_nodes))
                            child_nodes_pc = [self.ftn_pc[ftn_i] for ftn_i in child_nodes if ftn_i in self.ftn_pc.keys()]
                            self.fdg_pc[ftn_idx]=sorted(child_nodes_pc)

                # remove states, the functions of which have no dependent functions.
                if len(prt_no_child)>0 and len(prt_no_child)<len(prt_ftns) and len(prt_ftns)>=2:# have at least one function not lead node
                    for state in laserEVM.open_states:
                        # if not state.constraints.is_possible: continue
                        if state.node.function_name in self.FDG.ftn_to_index.keys():
                            ftn_idx = self.FDG.ftn_to_index[state.node.function_name]
                            if ftn_idx in prt_no_child:
                                laserEVM.open_states.remove(state)

            # prepare data for level-1 sequence execution
            if self._depth_==fdg.FDG_global.depth_all_ftns_reached+1:
                self._level_control_list.append(self._level_ + 1)
                self._level_+=1
                self._level_control_index += 1

                ftn_reached = self.OS_states[self._depth_-1].keys()
                ftn_uncovered = [ftn_idx for ftn_idx, item in enumerate(self.ftn_covered_mark) if item == 0]

                self.info_1st[1] = {}
                if len(ftn_reached) > 0:
                    temp_pc = [self.ftn_pc[ftn_i] for ftn_i in ftn_uncovered if ftn_i in self.ftn_pc.keys()]
                    temp_pc.sort()
                    for idx in ftn_reached:
                        # on all possible states, execute uncovered functions
                        self.info_1st[1][idx] = temp_pc


            # generate Sequence object
            if self._depth_ == fdg.FDG_global.depth_all_ftns_reached + 2:
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

            if self._depth_>fdg.FDG_global.depth_all_ftns_reached+1:
                self._level_control_index+=1

                if len(self._level_control_list)==self._level_control_index:
                    print(f'end of all levels')
                    fdg.FDG_global.transaction_count=self._depth_
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
                    if len(info_1st_temp)>0:
                        states_for_seq_execution =[]
                        for depth,ftn_dict in info_1st_temp.items():
                            for key_ftn in ftn_dict.keys():
                                # the states at depth_all_ftns_reached are already there.
                                if depth in self.OS_states.keys():
                                    if key_ftn in self.OS_states[depth].keys():
                                        states_for_seq_execution += self.OS_states[depth][key_ftn]


                        # place states to laserEVM instance
                        if len(states_for_seq_execution)>0:
                            laserEVM.open_states = states_for_seq_execution

                            # replace functions to be executd with its associated PCs
                            info_1={}
                            for depth,values in info_1st_temp.items():
                                for k,ftn_list in values.items():
                                    pc_values=[]
                                    for ftn_idx in ftn_list:
                                        if ftn_idx in self.ftn_pc.keys():
                                            pc_values.append(self.ftn_pc[ftn_idx])
                                    pc_values.sort()
                                    if k in info_1.keys():
                                        pc_list=info_1[k]+pc_values
                                        pc_list.sort()
                                        info_1[k]=pc_list
                                    else:info_1[k]=pc_values
                            self.info_1st[self._level_]=info_1

                            if info_following_temp:
                                for step,values in info_following_temp.items():
                                    self._level_control_list.append(-step) #  indicate the case to continue sequence execution
                                    for k,ftn_list in values.items():
                                        pc_values=[]
                                        for ftn_idx in ftn_list:
                                            if ftn_idx in self.ftn_pc.keys():
                                                pc_values.append(self.ftn_pc[ftn_idx])
                                        pc_values.sort()
                                        info_following_temp[step][k] = pc_values

                                self.info_following[self._level_]=info_following_temp
                            # prepare for the next level value
                            self._level_control_list.append(self._level_ + 1)

                        else:
                            # go the next level
                            self._level_control_list.append(self._level_ + 1)
                            self.info_1st[self._level_] = {}
                            print(f'no states for sequence execution at level: {self._level_}')
                    else:

                        self.info_1st[self._level_] = {}
                        # prepare the next level
                        self._level_control_list.append(self._level_+1)

                else:
                    # do not do anything here, just continue sequence execution
                    pass


        @symbolic_vm.laser_hook("stop_sym_trans_laserEVM")
        def stop_sym_trans_hook_laserEVM(laserEVM:LaserEVM):
            """
            - save states at the end of symbolic transactions for sequence execution
            - only need to states from depth 1 to fdg.FDG_global.depth_all_ftns_reached -1
            - the states at fdg.FDG_global.depth_all_ftns_reached are already saved in laserEVM
            - some saved states are used as initial states in sequence execution

            :param laserEVM:
            :return:
            """
            # assume that states are annotated by function sequences from which they are generated

            # save states from depth 2 to self._depth__all_ftns_reached
            if self._depth_<=fdg.FDG_global.depth_all_ftns_reached+1 and self._depth_>=1:
                    for state in laserEVM.open_states:
                        if not state.constraints.is_possible:continue
                        if state.node.function_name in self.FDG.ftn_to_index.keys():
                            ftn_name=state.node.function_name
                            ftn_idx=self.FDG.ftn_to_index[ftn_name]

                            if self._depth_==1:
                                self.ftn_start_nodes.append(ftn_idx)

                            self.ftn_covered_mark[ftn_idx] = 1# record the function covered                          

                            if len(self.ftn_wo_edges_not_covered)>0:
                                if ftn_idx in self.ftn_wo_edges_not_covered:
                                    self.ftn_wo_edges_not_covered.remove(ftn_idx)

                            
                            if ftn_idx not in self.OS_states[self._depth_].keys():
                                self.OS_states[self._depth_][ftn_idx] = [copy(state)]
                            else:
                                self.OS_states[self._depth_][ftn_idx] += [copy(state)]


            # record if the function covered
            if self._depth_ > fdg.FDG_global.depth_all_ftns_reached+1:
                for state in laserEVM.open_states:
                    if not state.constraints.is_possible: continue
                    if state.node.function_name in self.FDG.ftn_to_index.keys():
                        ftn_idx = self.FDG.ftn_to_index[state.node.function_name]
                        self.ftn_covered_mark[ftn_idx] = 1
                        if len(self.ftn_wo_edges_not_covered) > 0:
                            if ftn_idx in self.ftn_wo_edges_not_covered:
                                self.ftn_wo_edges_not_covered.remove(ftn_idx)

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
                    fdg.FDG_global.transaction_count =2



                # extract valid pc interval
                if len(self.gt_pc)>0:
                    self.valid_pc_interval = utils.get_valid_pc_interval(self.gt_pc, self.pc_control_interval[
                        'pc_interval_end'])

                # build FDG
                self.FDG.build_fdg_3d_array(self.ftn_start_nodes)

                # flag the case when there is no public function except fallback
                if self.FDG.num_ftn==1:
                    self.no_public_ftn=True
                self.all_ftn_pc_list.sort()
                # make sure 2 transactions are issued
                if self.FDG.depth_all_ftns_reached ==1:
                    # if self.ftn_covered_mark.count(0)==0:
                    #     fdg.FDG_global.transaction_count = 1 # stop at depth 1
                    # else:
                    fdg.FDG_global.transaction_count=2 # continue to explore one more step


                # update the value fdg.FDG_global.depth_all_ftns_reached
                if self.FDG.depth_all_ftns_reached>=2:
                    fdg.FDG_global.depth_all_ftns_reached=self.FDG.depth_all_ftns_reached


                # map function index to its pc
                for ftn_i in range(1,self.FDG.num_ftn):
                    selector=self.FDG.index_to_selector[ftn_i]
                    if selector in self.selector_pc.keys():
                        self.ftn_pc[ftn_i]=self.selector_pc[selector]
                        self.pc_ftn[self.selector_pc[selector]]=ftn_i


                # get values for ftn_wo_edges_not_covered_pc
                self.ftn_wo_edges_not_covered=self.FDG.nodes_wo_edges


            # get functions not covered
            if self._depth_<=fdg.FDG_global.depth_all_ftns_reached:
                if self.ftn_covered_mark.count(0)>0:
                    ftn_to_pc=[self.ftn_pc[ftn_i] for ftn_i, item in enumerate(self.ftn_covered_mark) if (item==0) and (ftn_i in self.ftn_pc.keys())]
                    ftn_to_pc.sort()
                    self.ftn_not_covered_pc_dict[self._depth_]=ftn_to_pc

            if self._depth_>=fdg.FDG_global.depth_all_ftns_reached:
                # stop symbolic execution if all functions are covered
                if self.ftn_covered_mark.count(0)==0:
                    # set to the current iteration, so that execution engine can stop
                    fdg.FDG_global.transaction_count = self._depth_


        #-------------------------------------------------
        ''' 
        - basically, this is done by:
            saving PCs for evaluating functions
            and changing machine state PC to PCs associated to specified functions
        '''
        #-------------------------------------------------

        @symbolic_vm.post_hook("DUP1")
        def dup1_hook(state: GlobalState):
            # if state.environment.active_function_name == 'fallback'\
            # and self._depth_ >= 2:
            if self._depth_ >= 2:
                # only consider DUP1 within a specified range
                pc_here = state.mstate.pc
                if len(self.gt_pc)==0:
                    if pc_here < self.pc_control_interval['pc_interval_start']: return
                    if pc_here > self.pc_control_interval['pc_interval_end']: return
                else:
                    if not utils.pc_is_valid(pc_here,self.valid_pc_interval):
                        return
                # handle the case when there is no public function except fallback
                if self.no_public_ftn:
                    pc_list = self.ftn_not_covered_pc_dict[self._depth_ - 1]
                    state.mstate.pc = utils.assign_pc_seq_exe_phase_dup1(
                        state.mstate.pc, self.all_ftn_pc_list,self.pc_control_interval)
                    return

                annotations = get_dependency_annotation(state)
                ftn_seq=annotations.ftn_seq
                if len(ftn_seq)==0: return
                pre_ftn_name=ftn_seq[-1] # get the function name

                # for FDG exploration
                if self._depth_ <= fdg.FDG_global.depth_all_ftns_reached:
                    state.mstate.pc=utils.assign_pc_fdg_phase_dup1(state.mstate.pc,
                                                                   pre_ftn_name,
                                                                   self.FDG.ftn_to_index,
                                                                   self.fdg_pc,
                                                                   self._depth_-1,
                                                                   self.ftn_not_covered_pc_dict,
                                                                   self.pc_control_interval)
                    return


                else: # for sequence execution
                    if pre_ftn_name not in self.FDG.ftn_to_index.keys():
                        state.mstate.pc=self.pc_control_interval['pc_interval_end']
                        return
                    pre_ftn_name_idx = self.FDG.ftn_to_index[pre_ftn_name]

                    # level-1 sequence execution
                    if self._depth_==fdg.FDG_global.depth_all_ftns_reached+1:
                        # pc_list=self.info_1st[1][pre_ftn_name_idx]
                        pc_list = self.ftn_not_covered_pc_dict[self._depth_-1]
                        state.mstate.pc = utils.assign_pc_seq_exe_phase_dup1(state.mstate.pc, pc_list,
                                                                             self.pc_control_interval)
                        return

                    # sequence execution with level>1
                    # for the first step
                    if self._level_==self._level_control_list[self._level_control_index]:
                        if pre_ftn_name_idx in self.info_1st[self._level_].keys():
                            pc_list = self.info_1st[self._level_][pre_ftn_name_idx]
                            state.mstate.pc = utils.assign_pc_seq_exe_phase_dup1(state.mstate.pc, pc_list,
                                                                                 self.pc_control_interval)
                            return
                        state.mstate.pc = self.pc_control_interval['pc_interval_end']
                        return

                    else: # for the following steps
                        # get the step for sequence execution
                        seq_execution_step=-self._level_control_list[self._level_control_index]
                        if seq_execution_step in self.info_following[self._level_].keys():
                            if pre_ftn_name_idx in self.info_following[self._level_][seq_execution_step].keys():
                                pc_list = self.info_following[self._level_][seq_execution_step][pre_ftn_name_idx]
                                state.mstate.pc = utils.assign_pc_seq_exe_phase_dup1(state.mstate.pc, pc_list,
                                                                                     self.pc_control_interval)
                                return
                        state.mstate.pc = self.pc_control_interval['pc_interval_end']
                        return

        @symbolic_vm.pre_hook("PUSH4")
        def push4_hook(state: GlobalState):
            # get the pc of opcode starting matching for each function
            # if state.environment.active_function_name == 'fallback' \
            if self._depth_ == 1:
                # assume that self.pc_control_interval is extracted
                if state.mstate.pc< self.pc_control_interval['pc_interval_start']:
                    return
                if 'pc_interval_end' in self.pc_control_interval.keys():
                    if state.mstate.pc > self.pc_control_interval['pc_interval_end']:
                        return
                ftn_selector = state.instruction['argument'][2:]
                if ftn_selector not in self.selector_pc.keys():
                    self.selector_pc[ftn_selector] = state.mstate.pc
                    if ftn_selector!='ffffffff':
                        self.all_ftn_pc_list.append(state.mstate.pc)


        @symbolic_vm.pre_hook("GT")
        def gt_hook(state: GlobalState):
            # get the pc of jumpdest, the start opcode of a block meaning the end of function mapping

            if self._depth_ == 1:
                if state.mstate.pc< self.pc_control_interval['pc_interval_start']:
                    return
                if 'pc_interval_end' in self.pc_control_interval.keys():
                    if state.mstate.pc > self.pc_control_interval['pc_interval_end']:
                        return

                self.gt_pc.append(state.mstate.pc)


        @symbolic_vm.pre_hook("CALLDATALOAD")
        def calldataload_hook(state: GlobalState):
            '''
            get the start pc for the valid pc interval
            '''
            if self._depth_==1:
                if 'pc_interval_start' not in self.pc_control_interval.keys():
                    self.pc_control_interval['pc_interval_start']=state.mstate.pc

        @symbolic_vm.pre_hook("CALLDATASIZE")
        def calldatasize_hook(state: GlobalState):
                if self._depth_ == 1:
                    if 'pc_signal_start' not in self.pc_control_interval.keys():
                        self.pc_control_interval['pc_signal_start'] = state.mstate.pc

        @symbolic_vm.pre_hook("JUMPDEST")
        def jumpdest_hook(state: GlobalState):
            '''
            get the maximum pc for the valid pc interval
            '''
            # get the pc of jumpdest, the start opcode of a block meaning the end of function mapping
            # if state.environment.active_function_name == 'fallback':

            # assume that the first occurance of jumpdest is the entry point for block of revert or fallback function

            if self._depth_ == 1 and 'pc_signal_start' in self.pc_control_interval.keys() and 'pc_interval_end' not in self.pc_control_interval.keys():
                self.pc_control_interval['pc_interval_end'] = state.mstate.pc


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
            # get function pairs from depth 1 to fdg.FDG_global.depth_all_ftns_reached
            if self._depth_>=2 and self._depth_<=fdg.FDG_global.depth_all_ftns_reached+1:
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
                self._depth_ = 0
                self._depth_=0
                return



from copy import copy

from virtualenv.config.convert import NoneType

from fdg import utils
from fdg.FDG_3d_array import FDG
from fdg.funtion_info import Function_info
from fdg.sequence import Sequence

from mythril.laser.ethereum.svm import LaserEVM
from mythril.laser.plugin.interface import LaserPlugin
from mythril.laser.plugin.builder import PluginBuilder
from mythril.laser.plugin.plugins.coverage import coverage_plugin
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


        self.solidity=''
        self.contract=''
        self.FDG=None       
        self.fdg_pc={} # save the pc of instruction for matching each function

        self.selector_pc={}
        self.ftn_pc={}
        self.pc_ftn={}
        self.pc_control_interval={}
        self.gt_pc = []
        self.valid_pc_interval = []


        self.OS_states={}
        self.ftn_pairs={} # save (parent, child)

        self.ftn_not_covered = []
        self.ftn_not_covered_pc = []
        self.ftn_sv_not_covered_pc = []
        self.ftn_wo_edges_not_covered = []

        self.ftn_instructions_coverage_info={}
        self.ftn_instructions_indices={}
        self.ftn_identifiers={}

        self.seq_1_parent_considered = {}
        self.seq_object = None
        self.seq_depth=-1
        self.seq_depth_max=0
        self.seq_current_package={} 
        self.seq_ftn_no_sequences=[]
        self.flag_ftn_no_sequence=False
        self.seq_ftn_no_sequence_count=0



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

            self.ftn_instructions_indices=fdg.FDG_global.ftns_instr_indices
            for ftn_full_name, identifier in fdg.FDG_global.method_identifiers.items():
                # remve '(...)' from function signature,
                # use pure name as key because self.ftn_instructions_indices uses only pure name as key
                self.ftn_identifiers[str(ftn_full_name).split('(')[0]]=identifier

            for ftn, ftn_instr_list in fdg.FDG_global.ftns_instr_indices.items():
                # if ftn=='constructor' or ftn=='fallback':continue
                if ftn == 'constructor': continue
                self.ftn_instructions_coverage_info[ftn]=[0/len(ftn_instr_list),ftn_instr_list]


        @symbolic_vm.laser_hook("stop_sym_exec")
        def stop_sym_exec_hook():
            print(f'End of symbolic execution')
            for ftn, ftn_cov in self.ftn_instructions_coverage_info.items():
                print("{:.2f}% coverage for {}".format(ftn_cov[0], ftn))



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
            ...
            add states to laserEVM.open_states so that they can be used
            as base states in the next iteration of symbolic transaction
            :param laserEVM: instance of LaserEVM
            :return:
            """
            self._iteration_+=1
            # define variables to save states and function pairs
            if self._iteration_ <= fdg.FDG_global.depth_all_ftns_reached+1:
                self.OS_states[self._iteration_]={}
                if self._iteration_>=2:
                    self.ftn_pairs[self._iteration_]=[]

        @symbolic_vm.laser_hook("stop_sym_trans_laserEVM")
        def stop_sym_trans_hook_laserEVM(laserEVM:LaserEVM):
            """
            - save states
            - only need to states from depth 1 to fdg.FDG_global.depth_all_ftns_reached+1
            - the states at fdg.FDG_global.depth_all_ftns_reached are already saved in laserEVM
            - some saved states are used as initial states in sequence execution

            :param laserEVM:
            :return:
            """
            # assume that states are annotated by function sequences from which they are generated
            if len(laserEVM.open_states)==0:
                print(f'no states are generated at depth {self._iteration_}')

            # save states and check which functions are covered
            if self._iteration_<=fdg.FDG_global.depth_all_ftns_reached+1:
                #save states from depth 1 to self._iteration__all_ftns_reached+1
                for state in laserEVM.open_states:
                    if not state.constraints.is_possible:continue
                    if state.node.function_name in self.FDG.ftn_to_index.keys():
                        ftn_name=state.node.function_name
                        ftn_idx=self.FDG.ftn_to_index[ftn_name]
                        if ftn_idx not in self.OS_states[self._iteration_].keys():
                            self.OS_states[self._iteration_][ftn_idx] = [copy(state)]
                        else:
                            self.OS_states[self._iteration_][ftn_idx] += [copy(state)]
                    else:
                        # in case that funtion pure name is the same, but function signature is not the same
                        ftn_pure_name=str(state.node.function_name).split('(')[0]
                        if ftn_pure_name in self.FDG.ftn_0_to_index.keys():
                            ftn_idx=self.FDG.ftn_0_to_index[ftn_pure_name]
                            if ftn_idx not in self.OS_states[self._iteration_].keys():
                                self.OS_states[self._iteration_][ftn_idx] = [copy(state)]
                            else:
                                self.OS_states[self._iteration_][ftn_idx] += [copy(state)]



            # check the code coverage for each function
            instr_cov_record_list = fdg.FDG_global.ftns_instr_cov[1]
            if len(instr_cov_record_list)>0:
                instr_array = np.array(instr_cov_record_list)
                self.ftn_not_covered = []
                self.ftn_sv_not_covered_pc=[]
                for ftn, ftn_instr_cov in self.ftn_instructions_coverage_info.items():
                    if ftn_instr_cov[0] == 100: continue
                    status = instr_array[fdg.FDG_global.ftns_instr_indices[ftn]]
                    cov_instr = sum(status)
                    cov = cov_instr / float(len(status)) * 100
                    self.ftn_instructions_coverage_info[ftn] = [cov, status]
                    if cov < 98:
                        # functions of public state variables can not be in FDG
                        if ftn in self.FDG.ftn_0_to_index.keys():
                            self.ftn_not_covered.append(self.FDG.ftn_0_to_index[ftn])
                        else:
                            # public functions of state variables, not in FDG.
                            identifier=self.ftn_identifiers[ftn] if ftn in self.ftn_identifiers.keys() else '0000000'
                            ftn_pc=self.selector_pc[identifier] if identifier in self.selector_pc.keys() else self.pc_control_interval['pc_interval_end']
                            self.ftn_sv_not_covered_pc.append(ftn_pc)


            # build FDG after the first symbolic transaction
            # at depth 1, build FDG
            if self._iteration_ == 1:
                ftn_start_nodes = self.OS_states[self._iteration_].keys()
                if len(ftn_start_nodes) == 0:
                    fdg.FDG_global.transaction_count = 2

                # extract valid pc interval
                if len(self.gt_pc) > 0:
                    self.valid_pc_interval = utils.get_valid_pc_interval(self.gt_pc,
                                                                         self.pc_control_interval[
                                                                             'pc_interval_end'])

                # build FDG
                self.FDG.build_fdg_3d_array(ftn_start_nodes)



                # make sure 2 transactions are issued
                if self.FDG.depth_all_ftns_reached == 1:
                    if self.ftn_not_covered.count(0) == 0 and fdg.FDG_global.coverage[1] >= 98:
                        fdg.FDG_global.transaction_count = 1  # stop at depth 1
                    else:
                        fdg.FDG_global.transaction_count = 2  # continue to explore one more step

                # update the value fdg.FDG_global.depth_all_ftns_reached
                if self.FDG.depth_all_ftns_reached >= 2:
                    fdg.FDG_global.depth_all_ftns_reached = self.FDG.depth_all_ftns_reached

                # map function index to its pc
                for ftn_i in range(2, self.FDG.num_ftn):  # constructor and fallback do not have selector
                    selector = self.FDG.index_to_selector[ftn_i]
                    if selector in self.selector_pc.keys():
                        self.ftn_pc[ftn_i] = self.selector_pc[selector]
                        self.pc_ftn[self.selector_pc[selector]] = ftn_i

                # get values for ftn_wo_edges_not_covered_pc
                self.ftn_wo_edges_not_covered = copy(self.FDG.nodes_wo_edges)
                if 1 in self.ftn_wo_edges_not_covered:  # remove fallback function
                    self.ftn_wo_edges_not_covered.remove(1)



            # get PCs of functions not covered
            ftn_to_pc = [self.ftn_pc[ftn_i] for ftn_i in self.ftn_not_covered if
                         ftn_i in self.ftn_pc.keys()]
            # add those uncovered public functions of state variables
            ftn_to_pc+=self.ftn_sv_not_covered_pc
            self.ftn_not_covered_pc = sorted(ftn_to_pc)

            # check and stop if all functions are covered
            if self._iteration_>=fdg.FDG_global.depth_all_ftns_reached:
                if len(self.ftn_not_covered)==0 and len(self.ftn_sv_not_covered_pc)==0:
                    # set to the current iteration, so that execution engine can stop
                    fdg.FDG_global.transaction_count = self._iteration_
                    return

            # # check if the functions without sequences is covered or not
            # if self._iteration_ > fdg.FDG_global.depth_all_ftns_reached + 1:
            #     temp_keep = []
            #     for ftn_item in self.seq_ftn_no_sequences:
            #         if ftn_item in self.ftn_not_covered:
            #             temp_keep.append(ftn_item)
            #     self.seq_ftn_no_sequences = temp_keep

            # update and/or add data in self.fdg_pc for the next iteration in FDG exploration
            if self._iteration_ <=fdg.FDG_global.depth_all_ftns_reached:
                # check if no-edge functions covered or not
                temp_keep = []
                for ftn_idx in self.ftn_wo_edges_not_covered:
                    if ftn_idx in self.ftn_not_covered:
                        temp_keep += [ftn_idx]
                self.ftn_wo_edges_not_covered = temp_keep

                prt_no_child = []  # record funtions who do not have dependent functions
                prt_ftns = self.OS_states[self._iteration_].keys()
                for ftn_idx in prt_ftns:
                    child_nodes = copy(self.ftn_wo_edges_not_covered)  # always consider those functions having no edges
                    if ftn_idx in self.FDG.graph.keys():
                        child_nodes += self.FDG.graph[ftn_idx]
                    if len(child_nodes) == 0:
                        prt_no_child.append(ftn_idx)
                        continue

                    # update or add data in self.fdg_pc
                    if self._iteration_ <fdg.FDG_global.depth_all_ftns_reached:
                        child_nodes = list(set(child_nodes))
                        child_nodes_pc = [self.ftn_pc[ftn_i] for ftn_i in child_nodes if ftn_i in self.ftn_pc.keys()]
                        self.fdg_pc[ftn_idx] = sorted(child_nodes_pc)
                    else:  # prepare for level-1 data
                        target_childs = list(set(self.ftn_not_covered).intersection(set(child_nodes)))
                        if len(target_childs) == 0:
                            prt_no_child.append(ftn_idx)
                        else:
                            self.seq_1_parent_considered[ftn_idx] = sorted(
                                [self.ftn_pc[ftn_i] for ftn_i in target_childs if ftn_i in self.ftn_pc.keys()])

                # remove states, the functions of which have no dependent functions.
                if len(prt_no_child) > 0 and len(prt_no_child) < len(prt_ftns) and len(
                    prt_ftns) >= 2:  # have at least one function not lead node
                    for state in laserEVM.open_states:
                        # if not state.constraints.is_possible: continue
                        if state.node.function_name in self.FDG.ftn_to_index.keys():
                            ftn_idx = self.FDG.ftn_to_index[state.node.function_name]
                            if ftn_idx in prt_no_child:
                                laserEVM.open_states.remove(state)


            # create Sequence object
            if self._iteration_ == fdg.FDG_global.depth_all_ftns_reached + 1:
                # sequence genration
                # build new fdg
                self.FDG.build_fdg_3d_array_new(self.ftn_pairs)

                # create an Sequence object
                seq_object = Sequence(self.FDG, self.ftn_not_covered, 5)
                # get ready for sequence generation(parpare some data for sequence generation)
                seq_object.prepare_sequence_generation()
                self.seq_ftn_no_sequences = seq_object.get_all_sequences()
                self.seq_object = seq_object


            # sequence execution (when multiple parents are considered in sequence generation)
            if self._iteration_ >= fdg.FDG_global.depth_all_ftns_reached + 1:
                if self.seq_depth == self.seq_depth_max:
                    self.seq_depth = -1

                if self.seq_depth == -1:
                    while (self.seq_depth == -1):
                        self.seq_depth = 0
                        self.seq_depth_max = 0
                        self.seq_current_package = self.seq_object.get_package_dict(self.ftn_not_covered)

                        # no sequence is available
                        if len(self.seq_current_package.keys()) == 0 or len(
                            self.seq_current_package['locate_state']) == 0:
                            print(f'no sequence generated or available')
                            # stop
                            # flag to start handling functions without sequences generated
                            self.flag_ftn_no_sequence=True
                            break

                        # locate states
                        states_for_seq_execution = []
                        for value in self.seq_current_package['locate_state']:
                            if value[0] in self.OS_states.keys():
                                if value[1] in self.OS_states[value[0]].keys():
                                    states_for_seq_execution += self.OS_states[value[0]][value[1]]

                        if len(states_for_seq_execution) > 0:
                            laserEVM.open_states = states_for_seq_execution
                            # map ftn to its pc
                            for depth, ftn_dict in self.seq_current_package.items():
                                if depth != 'locate_state':
                                    for ftn_i, values in ftn_dict.items():
                                        self.seq_current_package[depth][ftn_i]=sorted([self.ftn_pc[ftn] for ftn in values if ftn in self.ftn_pc.keys()])

                                    if self.seq_depth_max < depth:
                                        self.seq_depth_max = depth
                        else:
                            # go to the next package
                            self.seq_depth = -1
                else: # continue execution in the same package
                    self.seq_depth += 1


            # handle functions with no sequences generated
            if self.flag_ftn_no_sequence:
                # only execute one time
                self.seq_ftn_no_sequence_count+=1
                if self.seq_ftn_no_sequence_count>1:
                    fdg.FDG_global.transaction_count = self._iteration_
                    return
                if len(self.seq_ftn_no_sequences)>0:
                    # load states from FDG-guided phase
                    states_to_load=[]
                    states=self.OS_states[fdg.FDG_global.depth_all_ftns_reached]
                    for values in states.values():
                        for state in values:
                            states_to_load.append(state)
                    laserEVM.open_states=states_to_load

                    self.seq_ftn_no_sequences_pc=sorted([self.ftn_pc[ftn_i] for ftn_i in self.seq_ftn_no_sequences if ftn_i in self.ftn_pc.keys()])

                    # allow one more iteration
                    fdg.FDG_global.transaction_count=self._iteration_+1

                else:
                    # end executiong
                    fdg.FDG_global.transaction_count = self._iteration_



        ''' 
             changing machine state PC to PCs associated to specified functions at depth >1
        '''
        #-------------------------------------------------

        @symbolic_vm.post_hook("DUP1")
        def dup1_hook(state: GlobalState):
            if self._iteration_ >= 2:
                # only consider DUP1 within a specified range
                pc_here = state.mstate.pc
                if len(self.selector_pc)==0:
                    return
                if len(self.gt_pc)==0:
                    if pc_here < self.pc_control_interval['pc_interval_start']: return
                    if pc_here > self.pc_control_interval['pc_interval_end']: return
                else:
                    if not utils.pc_is_valid(pc_here,self.valid_pc_interval):
                        return

                # get the index of function, execution of which generates the state used as initial state in this transaction
                annotations = get_dependency_annotation(state)
                ftn_seq=annotations.ftn_seq
                if len(ftn_seq)==0: return
                pre_ftn_name=ftn_seq[-1] # get the function name
                pre_ftn_name_idx=-1
                if pre_ftn_name not in self.FDG.ftn_to_index.keys():
                    if str(pre_ftn_name).__contains__('('):
                        if str(pre_ftn_name).split('(')[0] in self.FDG.ftn_0_to_index.keys():
                            pre_ftn_name_idx= self.FDG.ftn_0_to_index[str(pre_ftn_name).split('(')[0]]
                else:
                    pre_ftn_name_idx = self.FDG.ftn_to_index[pre_ftn_name]


                # prepare a list of PCs of functions to be assigned
                pc_list = []
                if pre_ftn_name_idx==-1:  # functions do not in FDG, but change states
                    print(f'function {pre_ftn_name} does not in FDG, but changes states at iteration {self._iteration_-1}')
                elif self.flag_ftn_no_sequence:
                    pc_list=self.seq_ftn_no_sequences_pc
                else:
                    if self._iteration_ <= fdg.FDG_global.depth_all_ftns_reached:
                        pc_list = self.fdg_pc[pre_ftn_name_idx] if pre_ftn_name_idx in self.fdg_pc.keys() else []

                    elif self._iteration_ == fdg.FDG_global.depth_all_ftns_reached+1:
                        pc_list=pc_list = self.seq_1_parent_considered[pre_ftn_name_idx] if pre_ftn_name_idx in self.seq_1_parent_considered.keys() else []
                    else:
                        if self.seq_depth in self.seq_current_package.keys() and \
                            pre_ftn_name_idx in self.seq_current_package[self.seq_depth].keys():
                            pc_list = self.seq_current_package[self.seq_depth][pre_ftn_name_idx]


                # assign pc
                state.mstate.pc = utils.assign_pc_seq_exe_phase_dup1(state.mstate.pc, pc_list,
                                                                     self.ftn_not_covered_pc,
                                                                     self.pc_control_interval)






        #-------------------------------------------------
        ''' 
             save PC for each callable function at depth 1
        '''
        #-------------------------------------------------

        @symbolic_vm.pre_hook("PUSH4")
        def push4_hook(state: GlobalState):
            if self._iteration_ == 1:
                # assume that self.pc_control_interval is extracted
                if state.mstate.pc< self.pc_control_interval['pc_interval_start']:
                    return
                if 'pc_interval_end' in self.pc_control_interval.keys():
                    if state.mstate.pc > self.pc_control_interval['pc_interval_end']:
                        return
                ftn_selector = state.instruction['argument'][2:]
                if ftn_selector not in self.selector_pc.keys():
                    self.selector_pc[ftn_selector] = state.mstate.pc


        # -------------------------------------------------
        ''' 
             capture the correct PC interval at depth 1 that dispatcher matches functions
        '''
        # -------------------------------------------------

        @symbolic_vm.pre_hook("GT")
        def gt_hook(state: GlobalState):
            # get the pc of jumpdest, the start opcode of a block meaning the end of function mapping

            if self._iteration_ == 1:
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
            if self._iteration_==1:
                if 'pc_interval_start' not in self.pc_control_interval.keys():
                    self.pc_control_interval['pc_interval_start']=state.mstate.pc

        @symbolic_vm.pre_hook("CALLDATASIZE")
        def calldatasize_hook(state: GlobalState):
                if self._iteration_ == 1:
                    if 'pc_signal_start' not in self.pc_control_interval.keys():
                        self.pc_control_interval['pc_signal_start'] = state.mstate.pc

        @symbolic_vm.pre_hook("JUMPDEST")
        def jumpdest_hook(state: GlobalState):
            '''
            get the maximum pc for the valid pc interval
            '''
            # pc should larger than self.pc_control_interval['pc_interval_start']
            # assume that the first occurance of jumpdest is the entry point for block of revert or fallback function

            if self._iteration_ == 1 and 'pc_signal_start' in self.pc_control_interval.keys() and 'pc_interval_end' not in self.pc_control_interval.keys():
                self.pc_control_interval['pc_interval_end'] = state.mstate.pc

        # -------------------------------------------------
        ''' 
             check states at the end of transactions to see to which function it belongs
        '''
        # -------------------------------------------------

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
            if self._iteration_>=2 and self._iteration_<=fdg.FDG_global.depth_all_ftns_reached+1:
                ftn_seq=get_dependency_annotation(state).ftn_seq
                ftn_idx_seq=[]
                for ftn_full_name in ftn_seq[-2:]:
                    if ftn_full_name in self.FDG.ftn_to_index.keys():
                        ftn_idx_seq.append(self.FDG.ftn_to_index[ftn_full_name])
                    else:
                        if str(ftn_full_name).count('('):
                            ftn_pure_name=str(ftn_full_name).split('(')[0]
                            if ftn_pure_name in self.FDG.ftn_0_to_index.keys():
                                ftn_idx_seq.append(self.FDG.ftn_0_to_index[ftn_pure_name])
                if len(ftn_idx_seq)<2:
                    # no need to add function pairs
                    # ( quite possible one of function is function not in FDG like fallback function)
                    return
                if (ftn_idx_seq[0],ftn_idx_seq[1]) not in self.ftn_pairs[self._iteration_]:
                    self.ftn_pairs[self._iteration_]+=[(ftn_idx_seq[0],ftn_idx_seq[1])]

        @symbolic_vm.laser_hook("add_world_state")
        def world_state_filter_hook(state: GlobalState):
            if isinstance(state.current_transaction, ContractCreationTransaction):
                # Reset iteration variable
                self._iteration_ = 0
                return



from copy import copy

from virtualenv.config.convert import NoneType

from fdg import utils
# from fdg.FDG_3d_array import FDG
from fdg.FDG import FDG
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
        self.function_mark=[]
        self.fdg_pc={} # save the pc of instruction for matching each function

        self.all_pc_list={}

        self.selector_pc={}
        self.ftn_pc={}
        self.pc_ftn={}
        self.pc_control_interval={}
        self.gt_pc = []
        self.valid_pc_interval = []
        self.executed_ftn_pc={}

        self.OS_states={}


        self.uncovered_functions = []
        self.uncovered_ftn_pc_list = []

        self.ftn_sv_not_covered_pc = []

        self.leaf_nodes_wo_DD_edges=[]
        self.uncovered_leaf_nodes_wo_DD_edges = []
        self.nodes_unable_to_assign=[]

        # coverage related
        self.ftn_instructions_coverage_info={}
        self.ftn_instructions_indices={}
        self.ftn_identifiers={}
        # self.instruction_list=[]

        # sequence related
        self.seq_1_parent_considered = {}
        self.seq_object = None
        self.seq_depth=-1
        self.seq_depth_max=0
        self.seq_current_package={}
        self.sequences={}




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
            # build FDG
            function_info = Function_info(self.solidity, self.contract)
            self.FDG= FDG(function_info.functions_dict_slither())

            self.function_mark=[False]*self.FDG.num_ftn
            self.function_mark[0]=True
            self.function_mark[1]=True
            for ftn_idx in self.FDG.nodes_wo_DD_edges:
                self.function_mark[ftn_idx]=True
            self.uncovered_leaf_nodes_wo_DD_edges=self.FDG.nodes_wo_DD_edges




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

            # # check the code coverage for each function
            # instr_cov_record_list = fdg.FDG_global.ftns_instr_cov[1]
            # if len(instr_cov_record_list)>0:
            #     instr_array = np.array(instr_cov_record_list)
            #     for ftn, ftn_instr_cov in self.ftn_instructions_coverage_info.items():
            #         if ftn_instr_cov[0] == 100: continue
            #         status = instr_array[fdg.FDG_global.ftns_instr_indices[ftn]]
            #         opcodes=self.instruction_list[fdg.FDG_global.ftns_instr_indices[ftn]]
            #         opcode_idx_not_covered=list(np.invert(status))
            #         opcodes_not_covered=opcodes[opcode_idx_not_covered]
            #         print(f'{ftn},not covered: {opcodes_not_covered}')
            #         if ftn=='mint':
            #             print(f'mint: opcodes:{opcodes}')
            #         if ftn == 'transferFrom':
            #             print(f'transferFrom: opcodes:{opcodes}')
            # all_instr_idx_not_covered=list(np.invert(instr_array))
            # all_instr_not_covered=self.instruction_list[all_instr_idx_not_covered]
            # print(f'contract: not covered:{all_instr_not_covered}')

        @symbolic_vm.laser_hook("stop_sym_trans")
        def execute_stop_sym_trans_hook():
            # ----------------------------
            if self._iteration_ == 1:
                # extract valid pc interval
                if len(self.gt_pc) > 0:
                    self.valid_pc_interval = utils.get_valid_pc_interval(self.gt_pc,
                                                                         self.pc_control_interval[
                                                                             'pc_interval_end'])
                # map function index to its pc
                for ftn_i in range(2, self.FDG.num_ftn):  # constructor and fallback do not have selector
                    selector = self.FDG.index_to_selector[ftn_i]
                    if selector in self.selector_pc.keys():
                        self.ftn_pc[ftn_i] = self.selector_pc[selector]
                        self.pc_ftn[self.selector_pc[selector]] = ftn_i


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
                self.executed_ftn_pc[self._iteration_]={}


        @symbolic_vm.laser_hook("stop_sym_trans_laserEVM")
        def stop_sym_trans_hook_laserEVM(laserEVM:LaserEVM):
            """
            - save states
            - only need to states from depth 1 to fdg.FDG_global.depth_all_ftns_reached+1
            - some saved states are used as initial states in sequence execution

            :param laserEVM:
            :return:
            """
            # ----------------------------
            # save states
            if self._iteration_<=fdg.FDG_global.depth_all_ftns_reached+1:
                #save states from depth 1 to self._iteration__all_ftns_reached+1
                for state in laserEVM.open_states:
                    if not state.constraints.is_possible:continue
                    ftn_name=state.node.function_name
                    ftn_idx=0
                    if ftn_name in self.FDG.ftn_to_index.keys():
                        ftn_idx=self.FDG.ftn_to_index[ftn_name]
                    else:
                        # in case that funtion pure name is the same, but function signature is not the same
                        ftn_pure_name = str(ftn_name).split('(')[0]
                        if ftn_pure_name in self.FDG.ftn_0_to_index.keys():
                            ftn_idx = self.FDG.ftn_0_to_index[ftn_pure_name]
                    if ftn_idx>0:
                        if ftn_idx not in self.OS_states[self._iteration_].keys():
                            self.OS_states[self._iteration_][ftn_idx] = [copy(state)]
                        else:
                            self.OS_states[self._iteration_][ftn_idx] += [copy(state)]






            # ----------------------------
            # compute the depth (<=4)
            if not all(self.function_mark):
                if self._iteration_<=3:
                    reached_ftns = self.OS_states[self._iteration_].keys()
                    if self._iteration_==1:
                        for ftn_idx in reached_ftns:
                            self.function_mark[ftn_idx] = True
                    children=[self.FDG.graph[ftn_idx] for ftn_idx in reached_ftns if ftn_idx in self.FDG.graph.keys()]
                    children=[item for sublist in children for item in sublist]
                    if len(children)==0:
                        # assign all uncovered functions
                        self.function_mark = [True] * self.FDG.num_ftn
                        fdg.FDG_global.depth_all_ftns_reached = self._iteration_ + 1
                    else:
                        for child in children:
                            self.function_mark[child] = True
                        if all(self.function_mark):
                            fdg.FDG_global.depth_all_ftns_reached=self._iteration_+1
                        else:
                            if self._iteration_==3:
                                self.nodes_unable_to_assign=[idx for idx in range(1,self.FDG.num_ftn) if self.function_mark[idx]==False]
                                self.function_mark=[True]*self.FDG.num_ftn
                                fdg.FDG_global.depth_all_ftns_reached = self._iteration_ + 1
            else:
                if fdg.FDG_global.depth_all_ftns_reached==5:
                    fdg.FDG_global.depth_all_ftns_reached=self._iteration_+1



            # ----------------------------
            # check the code coverage for each function
            instr_cov_record_list = fdg.FDG_global.ftns_instr_cov
            if len(instr_cov_record_list)>0:
                instr_array = np.array(instr_cov_record_list)
                self.uncovered_functions = []
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
                            self.uncovered_functions.append(self.FDG.ftn_0_to_index[ftn])
                        else:
                            # public functions of state variables, not in FDG.
                            identifier=self.ftn_identifiers[ftn] if ftn in self.ftn_identifiers.keys() else '0000000'
                            ftn_pc=self.selector_pc[identifier] if identifier in self.selector_pc.keys() else self.pc_control_interval['pc_interval_end']
                            self.ftn_sv_not_covered_pc.append(ftn_pc)




            # ----------------------------
            # get PCs of functions not covered
            ftn_to_pc = [self.ftn_pc[ftn_i] for ftn_i in self.uncovered_functions if
                         ftn_i in self.ftn_pc.keys()]
            # add those uncovered public functions of state variables
            ftn_to_pc+=self.ftn_sv_not_covered_pc
            self.uncovered_ftn_pc_list = sorted(ftn_to_pc)

            #----------------------------
            # check and stop if all functions are covered
            if self._iteration_>=fdg.FDG_global.depth_all_ftns_reached:
                if len(self.uncovered_functions)==0 and len(self.ftn_sv_not_covered_pc)==0:
                    # set to the current iteration, so that execution engine can stop
                    fdg.FDG_global.transaction_count = self._iteration_
                    return

            # ----------------------------
            # update and/or add data in self.fdg_pc for the next iteration in FDG exploration
            if self._iteration_ <=fdg.FDG_global.depth_all_ftns_reached:
                prt_no_children=[]
                # check if no-edge functions covered or not
                self.uncovered_leaf_nodes_wo_DD_edges=\
                [ftn_idx for ftn_idx in self.uncovered_leaf_nodes_wo_DD_edges if ftn_idx in self.uncovered_functions]
                # check if unable_to_assign function is covered or not
                self.nodes_unable_to_assign=[ftn_idx for ftn_idx in self.nodes_unable_to_assign if ftn_idx in self.uncovered_functions]

                prt_ftns = self.OS_states[self._iteration_].keys()

                for ftn_idx in prt_ftns:
                    child_nodes = copy(self.uncovered_leaf_nodes_wo_DD_edges)  # always consider those functions having no edges
                    child_nodes+=copy(self.nodes_unable_to_assign)
                    if ftn_idx in self.FDG.graph.keys():
                        child_nodes += self.FDG.graph[ftn_idx]
                    if len(child_nodes) == 0:
                        if ftn_idx>1: # do not include fallback
                            prt_no_children.append(ftn_idx)
                        continue

                    # update or add data in self.fdg_pc
                    if self._iteration_ <fdg.FDG_global.depth_all_ftns_reached:
                        child_nodes = list(set(child_nodes))
                        child_nodes_pc = [self.ftn_pc[ftn_i] for ftn_i in child_nodes if ftn_i in self.ftn_pc.keys()]
                        self.fdg_pc[ftn_idx] = sorted(child_nodes_pc)
                    else:  # prepare data for sequence execution when only individual parents are considered
                        target_childs = list(set(self.uncovered_functions).intersection(set(child_nodes)))
                        if len(target_childs) == 0:
                            self.seq_1_parent_considered[ftn_idx]=[]
                        else:
                            self.seq_1_parent_considered[ftn_idx] = sorted(
                                [self.ftn_pc[ftn_i] for ftn_i in target_childs if ftn_i in self.ftn_pc.keys()])

                # remove states
                if self._iteration_ <= fdg.FDG_global.depth_all_ftns_reached and len(prt_no_children)>0:


                    children = [self.FDG.graph[ftn_idx] for ftn_idx in prt_ftns if
                                ftn_idx in self.FDG.graph.keys()]
                    children = [item for sublist in children for item in sublist]
                    children = list(set(children))
                    if len(children)>0:
                        for state in laserEVM.open_states:
                            ftn_name = state.node.function_name
                            ftn_idx = 0
                            if ftn_name in self.FDG.ftn_to_index.keys():
                                ftn_idx = self.FDG.ftn_to_index[ftn_name]
                            else:
                                # in case that funtion pure name is the same, but function signature is not the same
                                ftn_pure_name = str(ftn_name).split('(')[0]
                                if ftn_pure_name in self.FDG.ftn_0_to_index.keys():
                                    ftn_idx = self.FDG.ftn_0_to_index[ftn_pure_name]
                            if ftn_idx in prt_no_children:
                                laserEVM.open_states.remove(state)



            # create Sequence object
            if self._iteration_ == fdg.FDG_global.depth_all_ftns_reached + 1:
                # create an Sequence object
                seq_object = Sequence(self.FDG, self.uncovered_functions,self.sequences, 5)

                have_sequences,self.seq_ftn_no_sequences = seq_object.generate_sequences()
                self.seq_object = seq_object
                if not have_sequences:
                    fdg.FDG_global.transaction_count=self._iteration_
                    return


            # sequence execution (when multiple parents are considered in sequence generation)
            if self._iteration_ >= fdg.FDG_global.depth_all_ftns_reached + 1:
                if self.seq_depth == self.seq_depth_max:
                    self.seq_depth = -1

                if self.seq_depth == -1:
                    while (self.seq_depth == -1):
                        self.seq_depth = 0
                        self.seq_depth_max = 0
                        self.seq_current_package = self.seq_object.get_package_dict(self.uncovered_functions)

                        # no sequence is available
                        if len(self.seq_current_package.keys()) == 0 or len(
                            self.seq_current_package['locate_state']) == 0:
                            fdg.FDG_global.transaction_count = self._iteration_
                            print(f'no seqequnece available')
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

            # if no states generated in sequence execution phase, go to the next package
            if self._iteration_ >fdg.FDG_global.depth_all_ftns_reached + 1:
                if len(laserEVM.open_states)==0:
                    # go the next sequence package
                    self.seq_depth=-1


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

                else:
                    if self._iteration_ <= fdg.FDG_global.depth_all_ftns_reached:
                        pc_list = self.fdg_pc[pre_ftn_name_idx] if pre_ftn_name_idx in self.fdg_pc.keys() else []

                    elif self._iteration_ == fdg.FDG_global.depth_all_ftns_reached+1:
                        pc_list=pc_list = self.seq_1_parent_considered[pre_ftn_name_idx] if pre_ftn_name_idx in self.seq_1_parent_considered.keys() else []
                    else:
                        if self.seq_depth in self.seq_current_package.keys() and \
                            pre_ftn_name_idx in self.seq_current_package[self.seq_depth].keys():
                            pc_list = self.seq_current_package[self.seq_depth][pre_ftn_name_idx]

                    executed_pc=[]
                    if pre_ftn_name_idx in self.executed_ftn_pc[self._iteration_].keys():
                        executed_pc=self.executed_ftn_pc[self._iteration_-1][pre_ftn_name_idx]

                    # assign pc
                state.mstate.pc = utils.assign_pc(state.mstate.pc, pc_list,executed_pc,
                                                                     self.uncovered_ftn_pc_list,
                                                                     self.pc_control_interval)






        #-------------------------------------------------
        ''' 
             save PC for each callable function at depth 1
        '''
        #-------------------------------------------------

        @symbolic_vm.pre_hook("PUSH4")
        def push4_hook(state: GlobalState):
            if self._iteration_ == 1:
                # if len(self.instruction_list) == 0:
                #     self.instruction_list = np.array(state.environment.code.instruction_list)

                # assume that self.pc_control_interval is extracted
                if 'pc_interval_start' not in self.pc_control_interval.keys():
                    return
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
            save function sequences
            :param state:
            """
            if self._iteration_ >= 1 and self._iteration_ <= fdg.FDG_global.depth_all_ftns_reached + 1:
                ftn_seq = get_dependency_annotation(state).ftn_seq
                ftn_idx_seq = []
                for ftn_full_name in ftn_seq:
                    if ftn_full_name in self.FDG.ftn_to_index.keys():
                        ftn_idx_seq.append(self.FDG.ftn_to_index[ftn_full_name])
                    else:
                        if str(ftn_full_name).count('('):
                            ftn_pure_name = str(ftn_full_name).split('(')[0]
                            if ftn_pure_name in self.FDG.ftn_0_to_index.keys():
                                ftn_idx_seq.append(self.FDG.ftn_0_to_index[ftn_pure_name])
                if len(ftn_idx_seq)>0:
                    if ftn_idx_seq[-1] in self.sequences.keys():
                        if ftn_idx_seq not in self.sequences[ftn_idx_seq[-1]]:
                            self.sequences[ftn_idx_seq[-1]]+=[ftn_idx_seq]
                    else:self.sequences[ftn_idx_seq[-1]]=[ftn_idx_seq]

                    if ftn_idx_seq[-1] not in self.executed_ftn_pc[self._iteration_].keys():
                        self.executed_ftn_pc[self._iteration_][ftn_idx_seq[-1]]=[self.ftn_pc[ftn_idx] for ftn_idx in ftn_idx_seq if ftn_idx in self.ftn_pc.keys()]




        @symbolic_vm.laser_hook("add_world_state")
        def world_state_filter_hook(state: GlobalState):
            if isinstance(state.current_transaction, ContractCreationTransaction):
                # Reset iteration variable
                self._iteration_ = 0
                return



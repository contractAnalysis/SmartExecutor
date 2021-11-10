# support FDG-guided execution and sequence execution

from copy import copy
from fdg import utils
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


class SSE_prunerBuilder(PluginBuilder):
    name = "sse-support sequence execution"

    def __call__(self, *args, **kwargs):
        return sse()


class sse(LaserPlugin):
    """ """

    def __init__(self):
        self._reset()

    def _reset(self):
        print(f'sequence={fdg.FDG_global.sequences}')
        self._iteration_ = 0
        self.solidity = ''
        self.contract = ''
        self.FDG = None


        self.function_mark = []  # used to record if a function is assigned or the execution of it succeeds at depth 1
        self.fdg_pc = {}  # assign chilren function entry PCs to each function used in FDG-guided execution phase

        # used to control the execution flow
        self.all_pc_list = {}
        self.selector_pc = {}
        self.ftn_pc = {}
        self.pc_ftn = {}
        self.pc_control_interval = {}
        self.gt_pc = []
        self.valid_pc_interval = []

        # save data during in FDG-guided execution phase
        self.OS_states = {}  # save in-between transaction states
        self.executed_ftn_pc = {}  # save executed ftn_pc sequences
        self.sequences = {}  # save vailid sequencess

        self.uncovered_functions = []
        self.uncovered_functions_pc_list = []
        self.ftn_special_pc = []
        self.uncovered_leaf_nodes_wo_DD_edges = []
        self.ftn_unable_to_assign = []

        # coverage related
        self.ftn_instructions_coverage_info = {}
        self.ftn_instructions_indices = {}
        self.ftn_identifiers = {}
        self.instruction_list = []

        # sequence related
        self.seq_object = None

        self.flag_no_sequence_generated_handle = False
        self.ftn_no_sequences_pc_list = []
        self.ftn_special_pc = []
        self.ftn_special_pc__no_sequence_pc_combined = []


        self.states_available = []
        self.states_available_depth = 0
        self.states_available_depth_index = 0


        self.sequences_given=[]
        self.sequences_given_pc=[]
        self.sequences_given_ftn_idx=[]
        self.sequences_given_length=0
        self.sequences_given_index=0

        self.cur_sequence = []
        self.cur_sequence_depth = 0
        self.cur_sequence_depth_index = 0



    def initialize(self, symbolic_vm: LaserEVM) -> None:
        """Initializes the FDG_pruner
        :param symbolic_vm
        """
        self._reset()

        @symbolic_vm.laser_hook("start_sym_exec")
        def start_sym_exec_hook():
            # initialize FDG
            self.solidity = fdg.FDG_global.solidity_path
            self.contract = fdg.FDG_global.contract

            # build FDG
            function_info = Function_info(self.solidity, self.contract)
            self.FDG = FDG(function_info.functions_dict_slither())

            # get sequences
            sequences=fdg.FDG_global.sequences.split(";")
            for seq in sequences:
                seq_list=seq.split(",")
                if len(seq_list)>=2:
                    [[1,seq_list[0]]]+seq_list[1:]
                self.sequences_given.append([[1,seq_list[0]]]+seq_list[1:])
            self.cur_sequence_depth=len(self.sequences_given)



            self.ftn_instructions_indices = fdg.FDG_global.ftns_instr_indices
            for ftn_full_name, identifier in fdg.FDG_global.method_identifiers.items():
                # remove '(...)' from function signature,
                # use pure name as key because self.ftn_instructions_indices uses only pure name as key
                self.ftn_identifiers[str(ftn_full_name).split('(')[0]] = identifier

            for ftn, ftn_instr_list in fdg.FDG_global.ftns_instr_indices.items():
                # if ftn=='constructor' or ftn=='fallback':continue
                if ftn == 'constructor': continue
                self.ftn_instructions_coverage_info[ftn] = [0 / len(ftn_instr_list), ftn_instr_list]

        @symbolic_vm.laser_hook("stop_sym_exec")
        def stop_sym_exec_hook():
            # compute coverage
            instr_cov_record_list = fdg.FDG_global.ftns_instr_cov
            if len(instr_cov_record_list) > 0:
                instr_array = np.array(instr_cov_record_list)
                for ftn, ftn_instr_cov in self.ftn_instructions_coverage_info.items():
                    if ftn_instr_cov[0] == 100: continue
                    status = instr_array[fdg.FDG_global.ftns_instr_indices[ftn]]
                    cov_instr = sum(status)
                    cov = cov_instr / float(len(status)) * 100
                    self.ftn_instructions_coverage_info[ftn] = [cov, status]


            if fdg.FDG_global.print_ftn_coverage == 1:
                print(f'End of symbolic execution')
                for ftn, ftn_cov in self.ftn_instructions_coverage_info.items():
                    print("{:.2f}% coverage for {}".format(ftn_cov[0], ftn))

            # # check the code coverage for each function
            # instr_cov_record_list = fdg.FDG_global.ftns_instr_cov
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

                # map fallback function to the max pc in pc_control_interval
                if 'pc_interval_end' in self.pc_control_interval.keys():
                    self.ftn_pc[1] = self.pc_control_interval['pc_interval_end']
                    self.pc_ftn[self.pc_control_interval['pc_interval_end']] = 1



        # -------------------------------------------------
        ''' 
          new hook methods for changing laserEVM instance
        - add states 
        - save states        
        '''

        # -------------------------------------------------
        @symbolic_vm.laser_hook("start_sym_trans_laserEVM")
        def start_sym_trans_hook_laserEVM(laserEVM: LaserEVM):
            """
            ...
            add states to laserEVM.open_states so that they can be used
            as base states in the next iteration of symbolic transaction
            :param laserEVM: instance of LaserEVM
            :return:
            """
            self._iteration_ += 1
            if self._iteration_==1:
                self.OS_states[self._iteration_]={}

            if self._iteration_==2:
                # initialize sequence execution
                self.cur_sequence = self.sequences_given_ftn_idx[0]
                self.cur_sequence_depth = len(self.cur_sequence)
                if self.cur_sequence_depth >= 2:
                    self.cur_sequence_depth_index = 1
                    laserEVM.open_states = copy(self.OS_states[1][self.cur_sequence[0][1]])
                    print(
                        f'iteration {self._iteration_},\texecute function: {self.sequences_given[self.sequences_given_index][self.cur_sequence_depth_index]}')

                else:
                    laserEVM.open_states = []
                    self.cur_sequence_depth_index = self.cur_sequence_depth
            elif self._iteration_>2:
                if self.cur_sequence_depth_index == self.cur_sequence_depth-1:

                    if self.sequences_given_index==self.sequences_given_length-1:
                        fdg.FDG_global.transaction_count=self._iteration_
                        laserEVM.open_states=[]

                    else:
                        # execute the next sequence
                        if self.sequences_given_index<self.sequences_given_length-1:
                            self.sequences_given_index+=1
                            self.cur_sequence=self.sequences_given_ftn_idx[self.sequences_given_index]
                            self.cur_sequence_depth=len(self.cur_sequence)
                            self.cur_sequence_depth_index = 1
                            laserEVM.open_states=copy(self.OS_states[1][self.cur_sequence[0][1]])
                        else:
                            fdg.FDG_global.transaction_count = self._iteration_
                            laserEVM.open_states = []

                else:
                    self.cur_sequence_depth_index+=1
                    print(f'iteration {self._iteration_},\texecute function: {self.sequences_given[self.sequences_given_index][self.cur_sequence_depth_index]}')







        @symbolic_vm.laser_hook("stop_sym_trans_laserEVM")
        def stop_sym_trans_hook_laserEVM(laserEVM: LaserEVM):
            """
            - save states
            - only need to states from depth 1 to fdg.FDG_global.depth_all_ftns_reached+1
            - some saved states are used as initial states in sequence execution

            :param laserEVM:
            :return:
            """

            if self._iteration_==1:
                # ----------------------------
                # save states

                for state in laserEVM.open_states:
                    if not state.constraints.is_possible: continue
                    ftn_name = state.node.function_name
                    ftn_idx = 0
                    if ftn_name in self.FDG.ftn_to_index.keys():
                        ftn_idx = self.FDG.ftn_to_index[ftn_name]
                    else:
                        # in case that funtion pure name is the same, but function signature is not the same
                        ftn_pure_name = str(ftn_name).split('(')[0]
                        if ftn_pure_name in self.FDG.ftn_0_to_index.keys():
                            ftn_idx = self.FDG.ftn_0_to_index[ftn_pure_name]
                    if ftn_idx > 0:
                        # save states at FDG-guided execution phase
                        if ftn_idx not in self.OS_states[self._iteration_].keys():
                            self.OS_states[self._iteration_][ftn_idx] = [copy(state)]
                        else:
                            self.OS_states[self._iteration_][ftn_idx] += [copy(state)]

                #-------------------------------
                # prepare for the following steps: save function indices and corresponding PCs
                for seq in self.sequences_given:
                    seq_pc = []
                    seq_idx=[]
                    for ftn_no_parameter in seq:
                        ftn_name = ftn_no_parameter
                        flag = False
                        if len(ftn_no_parameter) == 2:
                            ftn_name = ftn_no_parameter[1]
                            flag = True

                        if ftn_name in self.FDG.ftn_0_to_index.keys():
                            ftn_idx = self.FDG.ftn_0_to_index[ftn_name]
                            pc = self.ftn_pc[ftn_idx] if ftn_idx in self.ftn_pc.keys() else self.pc_control_interval[
                                "pc_interval_end"]
                            if flag:
                                seq_pc.append([1, pc])
                                seq_idx.append([1,ftn_idx])
                                flag = False
                            else:
                                seq_pc.append(pc)
                                seq_idx.append(ftn_idx)

                        else:
                            print(f'{ftn_name} does not in FDG!')

                    self.sequences_given_pc.append(seq_pc)
                    self.sequences_given_ftn_idx.append(seq_idx)





        ''' 
             changing machine state PC to PCs associated to specified functions at depth >1
        '''

        # -------------------------------------------------

        @symbolic_vm.post_hook("DUP1")
        def dup1_hook(state: GlobalState):
            if self._iteration_ >= 2:

                # only consider DUP1 within a specified range
                pc_here = state.mstate.pc
                if len(self.selector_pc) == 0:
                    return
                if len(self.gt_pc) == 0:
                    if pc_here < self.pc_control_interval['pc_interval_start']: return
                    if pc_here > self.pc_control_interval['pc_interval_end']: return
                else:
                    if not utils.pc_is_valid(pc_here, self.valid_pc_interval):
                        return

                if self.sequences_given_index==self.sequences_given_length or \
                    self. cur_sequence_depth_index==self.cur_sequence_depth:
                    pc_list=[]
                else:
                    pc_list = [self.sequences_given_pc[self.sequences_given_index][self.cur_sequence_depth_index]]

                # assign pc
                state.mstate.pc = utils.assign_pc_fdg_phase_dup1(state.mstate.pc, pc_list,
                                                                 self.pc_control_interval)



        # -------------------------------------------------
        ''' 
             save PC for each callable function at depth 1
        '''

        # -------------------------------------------------

        @symbolic_vm.pre_hook("PUSH4")
        def push4_hook(state: GlobalState):
            if self._iteration_ == 1:
                # if len(self.instruction_list) == 0:
                #     self.instruction_list = np.array(state.environment.code.instruction_list)

                # assume that self.pc_control_interval is extracted
                if 'pc_interval_start' not in self.pc_control_interval.keys():
                    return
                if state.mstate.pc < self.pc_control_interval['pc_interval_start']:
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
                if state.mstate.pc < self.pc_control_interval['pc_interval_start']:
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
            if self._iteration_ == 1:
                if 'pc_interval_start' not in self.pc_control_interval.keys():
                    self.pc_control_interval['pc_interval_start'] = state.mstate.pc

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

            # get valid sequences from states
            if self._iteration_ >= 1 and self._iteration_ <= fdg.FDG_global.depth_all_ftns_reached \
                or self.cur_sequence_depth_index >= 1 and self.cur_sequence_depth_index == self.cur_sequence_depth - 1:

                seq = _get_valid_sequence_from_state(state)
                if len(seq) >= 1:
                    if seq[-1] in self.sequences.keys():
                        if seq not in self.sequences[seq[-1]]:
                            self.sequences[seq[-1]] += [seq]
                    else:
                        self.sequences[seq[-1]] = [seq]

        def _get_valid_sequence_from_state(state: GlobalState):
            """
            get valid sequences from global states
            :param state:
            """
            ftn_seq = get_dependency_annotation(state).ftn_seq
            ftn_idx_seq = []
            for ftn_full_name in ftn_seq:
                if ftn_full_name in self.FDG.ftn_to_index.keys():
                    ftn_idx_seq.append(self.FDG.ftn_to_index[ftn_full_name])
                else:
                    ftn_pure_name = ftn_full_name
                    if str(ftn_full_name).count('('):
                        ftn_pure_name = str(ftn_full_name).split('(')[0]
                    if ftn_pure_name in self.FDG.ftn_0_to_index.keys():
                        ftn_idx_seq.append(self.FDG.ftn_0_to_index[ftn_pure_name])
            return ftn_idx_seq

        def _print_state_info(state:GlobalState)->None:
            print(f'==== constraints ====')
            for constraint in state.world_state.constraints:
                print(f'\t {constraint}')
            print(f'==== state.environment.active_account ====')
            print(f'\t {state.environment.active_account.address}')

            print(f'==== storage of the active_account ====')
            for key, value in state.environment.active_account.storage.printable_storage.items():
                print(f'\t key {key}  value {value}')

            print(f'==== memory ====')
            mem_size = state.mstate.memory_size
            for i in range(int(mem_size / 32)):
                word = state.mstate.memory.get_word_at(i)
                print(f'\t {word}')
            print(f'==== stack ====')
            for item in state.mstate.stack:
                print(f'\t {item}')





        @symbolic_vm.laser_hook("add_world_state")
        def world_state_filter_hook(state: GlobalState):
            if isinstance(state.current_transaction, ContractCreationTransaction):
                # Reset iteration variable
                self._iteration_ = 0
                return


    def _request_next_sequence(self, laserEVM: LaserEVM):
        self.cur_sequence=self.seq_object.get_one_sequence(self.uncovered_functions)

        if len(self.cur_sequence) < 2:  # the generated sequence has length >=2
            self.flag_sequence_handle = False
            # get functions that no sequences are generated for
            self.ftn_no_sequences_pc_list = [self.ftn_pc[ftn_idx] for ftn_idx in
                                             self.seq_object.ftn_no_sequences if
                                             ftn_idx in self.ftn_pc.keys()]
            self.ftn_no_sequences_pc_list.sort()

            if len(self.ftn_no_sequences_pc_list) > 0 or len(self.ftn_special_pc) > 0:
                self.flag_no_sequence_generated_handle = True
                self.states_available_depth_index = self.states_available_depth = 0
                self.states_available = []
        else:
            self.cur_sequence_depth = len(self.cur_sequence)
            self.cur_sequence_depth_index = 1
            # get states for the sequence
            # get the states for this sequence
            depth = self.cur_sequence[0][0]
            ftn_idx = self.cur_sequence[0][1]
            if depth in self.OS_states.keys():
                if ftn_idx in self.OS_states[depth].keys():
                    self.states_available = self.OS_states[depth][ftn_idx]
            self.states_available_depth = len(self.states_available)
            laserEVM.open_states=[copy(self.states_available[self.states_available_depth_index])]
    def _update_coverage(self):
        instr_cov_record_list = fdg.FDG_global.ftns_instr_cov
        if len(instr_cov_record_list) > 0:
            instr_array = np.array(instr_cov_record_list)
            self.uncovered_functions = []
            self.ftn_special_pc = []
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
                        # FDG is empty
                        # or public functions of state variables, not in FDG.
                        # if len(self.FDG.nodes)==0: # we do not consider function of state variables
                        identifier = self.ftn_identifiers[ftn] if ftn in self.ftn_identifiers.keys() else '0000000'
                        ftn_pc = self.selector_pc[identifier] if identifier in self.selector_pc.keys() else \
                            self.pc_control_interval['pc_interval_end']
                        self.ftn_special_pc.append(ftn_pc)

from copy import copy
from itertools import permutations, combinations
import numpy as np
import itertools as it

import fdg
from fdg import lists_merge, utils
from fdg.FDG import FDG
from fdg.funtion_info import Function_info


class Sequence():
    def __init__(self, fdg=None, ftn_idx_not_covered=[], valid_sequences={}, num_seq_limit=5):
        self.fdg = fdg
        self.ftn_idx_not_covered = ftn_idx_not_covered  # target function
        self.all_valid_sequences = {}
        # initialize all_valid_sequences
        for key, value in valid_sequences.items():
            if len(value) > 0:
                ele_size = [len(item) for item in value]
                min_size = min(ele_size)
                shortest_seq = [item for item in value if len(item) == min_size]
                self.all_valid_sequences[key] = {'sequences': value, 'seq_num': len(value), 'shortest': shortest_seq,
                                                 'shortest_depth': min_size}
            else:
                self.all_valid_sequences[key] = {'sequences': [], 'seq_num': 0, 'shortest': [], 'shortest_depth': 0}

        self.num_seq_limit = num_seq_limit

        self.generated_sequences = {}
        self.generated_sequences_assignment = []
        self.ftn_idx_to_column_idx = {}
        self.column_idx_to_ftn_idx = {}
        self.package_list = []

        self.ftn_no_sequences = []
        self.flag_ftn_no_sequences = False

    def _get_label_parent(self, ftn_idx) -> dict:
        """
        get parents and their labels based on the write and read of state variables
        :param ftn_idx:
        :return:dict: key: label index, value: parent indices
        """
        l_p_dict = {}
        sv_read_ary = self.fdg.sv_read[ftn_idx, :]
        sv_read = sv_read_ary[sv_read_ary >= 0]
        for sv in sv_read:
            ftn_write_ary = self.fdg.sv_write[:, sv]
            ftn_write = np.where(ftn_write_ary >= 0)[0]
            l_p_dict[sv] = [ftn_i for ftn_i in ftn_write if ftn_i != ftn_idx]
        return l_p_dict

    # consider all sequences for each parent
    def _get_sequences_by_level_3_4(self, parent_groups: list, num_groups: int, ftn_idx):
        """
        consider all sequences
        :param parent_groups:[[2,3,5],[3,6]]
        :param num_groups : number of parents to consider
        :param ftn_idx : the child
        :return:
        """
        # assert self.ftn_seq_and_shortest_seq_dict

        # get sequences through parent combination
        parent_combinations = utils.get_combination(parent_groups, num_groups)

        collection_seq = []  # save the sequences generated
        for comb in parent_combinations:
            # check if each parent has shortest sequence or not

            p_seq_num = []
            for p_element in comb:
                if p_element in self.all_valid_sequences.keys():
                    p_seq_num.append(self.all_valid_sequences[p_element]['seq_num'])
                else:
                    p_seq_num.append(0)

            # if one parent does not have sequence, ignore this parent combination
            if p_seq_num.count(0) > 0: continue

            # some parent has multiple shortest sequences, so need to do combination
            p_seq_idx_list = [list(range(num)) for num in p_seq_num]
            p_seq_index_comb = [list(com) for com in it.product(*p_seq_idx_list)]
            for p_seq_index in p_seq_index_comb:
                # replace each parent with its shortest sequence
                # (remove 0: constructor is ignored.  reverse: so that parent itself is the last element in its shortest sequence )
                sequence_list = [
                    self.all_valid_sequences[p_ele]['sequences'][index] \
                    for p_ele, index in zip(comb, p_seq_index)]
                if fdg.FDG_global.control_level == 3:
                    merge_seq = self._merge_sequences_ordered(sequence_list)
                else:
                    merge_seq = self._merge_sequences_all_permutation(sequence_list)

                for seq in merge_seq:
                    if len(seq) == 1: continue
                    temp = seq + [ftn_idx]
                    if temp not in collection_seq:
                        collection_seq.append(temp)
        return collection_seq

    # consider only the shortest sequences for each parent
    def _get_sequences_by_level_1_2(self, parent_groups: list, num_groups: int, ftn_idx):
        """
        Only consider shortest sequences
        :param parent_groups:[[2,3,5],[3,6]]
        :param num_groups : number of parents to consider
        :param ftn_idx : the child
        :return:
        """

        # get sequences through parent combination
        parent_combinations = utils.get_combination(parent_groups, num_groups)

        collection_seq = []  # save the sequences generated
        for comb in parent_combinations:
            p_seq_num = []
            for p_element in comb:
                if p_element in self.all_valid_sequences.keys():
                    p_seq_num.append(len(self.all_valid_sequences[p_element]['shortest']))
                else:
                    p_seq_num.append(0)

            # if one parent does not have sequence, ignore this parent combination
            if p_seq_num.count(0) > 0: continue

            # some parent has multiple shortest sequences, so need to do combination
            p_seq_idx_list = [list(range(num)) for num in p_seq_num]
            p_seq_index_comb = [list(com) for com in it.product(*p_seq_idx_list)]
            for p_seq_index in p_seq_index_comb:
                # replace each parent with its shortest sequence
                # (remove 0: constructor is ignored.  reverse: so that parent itself is the last element in its shortest sequence )
                sequence_list = [
                    self.all_valid_sequences[p_ele]['shortest'][index] \
                    for p_ele, index in zip(comb, p_seq_index)]

                if fdg.FDG_global.control_level == 1:
                    merge_seq = self._merge_sequences_ordered(sequence_list)
                else:
                    merge_seq = self._merge_sequences_all_permutation(sequence_list)

                for seq in merge_seq:
                    if len(seq) == 1: continue
                    temp = seq + [ftn_idx]
                    if temp not in collection_seq:
                        collection_seq.append(temp)

        return collection_seq

    def _merge_sequences_all_permutation(self, nested_list: list):
        """
        each permulation of nested list is created one sequence
        the prefix of the first sequence is fixed
        :param nested_list:
        :return:
        """
        result = []
        ele_num = len(nested_list)
        if ele_num == 1: return []

        permulation_nested_list = permutations(nested_list)
        for per_tuple in permulation_nested_list:
            temp_list = list(per_tuple)
            p_first = temp_list[0]
            p_first_length = len(temp_list[0])
            merge_seq = p_first
            for i in range(1, ele_num):
                merge_seq = lists_merge.merge_fix_list_1_specified_lenth_no_repeated_elements(merge_seq, temp_list[i],
                                                                                              p_first_length)
            # form final sequence
            # remove the elements from the first parent
            final_seq = [[p_first_length, p_first[-1]]]
            final_seq += merge_seq[p_first_length:]
            if final_seq not in result:
                result.append(final_seq)
        return result

    def _merge_sequences_ordered(self, nested_list: list):
        """
        the prefix of the first sequence is fixed
        :param nested_list:
        :return:
        """

        ele_num = len(nested_list)
        if ele_num == 1: return []
        # order elements in nested_list based on len
        order_nested_list = sorted(nested_list, key=len, reverse=True)

        # merge sequence
        p_first = order_nested_list[0]
        p_first_length = len(order_nested_list[0])
        merge_seq = p_first

        for i in range(1, ele_num):
            merge_seq = lists_merge.merge_fix_list_1_specified_lenth_no_repeated_elements(merge_seq,
                                                                                          order_nested_list[i],
                                                                                          p_first_length)

        # form final sequence
        # remove the elements from the first parent
        final_seq = [[p_first_length, p_first[-1]]]
        final_seq += merge_seq[p_first_length:]
        return [final_seq]

        # consider only the shortest sequences for each parent

    def _get_sequences_by_level_0(self, parent_groups: list, num_groups: int, ftn_idx):
        """
        Only consider shortest sequences
        :param parent_groups:[[2,3,5],[3,6]]
        :param num_groups : number of parents to consider
        :param ftn_idx : the child
        :return:
        """

        # get sequences through parent combination
        parent_combinations = utils.get_combination(parent_groups, num_groups)

        collection_seq = []  # save the sequences generated
        for comb in parent_combinations:
            p_seq_num = []
            for p_element in comb:
                if p_element in self.all_valid_sequences.keys():
                    p_seq_num.append(len(self.all_valid_sequences[p_element]['shortest']))
                else:
                    p_seq_num.append(0)

            # if one parent does not have sequence, ignore this parent combination
            if p_seq_num.count(0) > 0: continue

            # some parent has multiple shortest sequences, so need to do combination
            p_seq_idx_list = [list(range(num)) for num in p_seq_num]
            p_seq_index_comb = [list(com) for com in it.product(*p_seq_idx_list)]
            for p_seq_index in p_seq_index_comb:
                # replace each parent with its shortest sequence
                # (remove 0: constructor is ignored.  reverse: so that parent itself is the last element in its shortest sequence )
                sequence_list = [
                    self.all_valid_sequences[p_ele]['shortest'][index] \
                    for p_ele, index in zip(comb, p_seq_index)]
                merge_seq = self._merge_sequences(sequence_list)

                for seq in merge_seq:
                    if len(seq) == 1: continue
                    temp = seq + [ftn_idx]
                    if temp not in collection_seq:
                        collection_seq.append(temp)
        return collection_seq

    def _get_sequences_by_level_5(self, parents: list, ftn_idx):
        """
        consider all parents at one time
        use the shortest sequences to represent parents
        the prefix of the first sequence is not fixed
        :param parents:[2,3,5]

        :param ftn_idx : the child
        :return:
        """
        collection_seq=[]
        parents_having_seq=[]
        p_seq_num = []  # save number of parents' shortest seqeunces
        for p_element in parents:
            if p_element in self.all_valid_sequences.keys():
                p_seq_num.append(len(self.all_valid_sequences[p_element]['shortest']))
                parents_having_seq.append(p_element)


        # some parent has multiple shortest sequences, so need to do combination
        p_seq_idx_list = [list(range(num)) for num in p_seq_num]
        p_seq_index_comb = [list(com) for com in it.product(*p_seq_idx_list)]
        for p_seq_index in p_seq_index_comb:
            # replace each parent with its shortest sequence
            # (remove 0: constructor is ignored.  reverse: so that parent itself is the last element in its shortest sequence )
            sequence_list = [
                self.all_valid_sequences[p_ele]['shortest'][index] \
                for p_ele, index in zip(parents_having_seq, p_seq_index)]

            permulation_sequence_list= permutations(sequence_list)
            for seq_list in permulation_sequence_list:

                merge_seq = self._merge_sequences(seq_list)

                for seq in merge_seq:
                    if len(seq) == 1: continue
                    temp = seq + [ftn_idx]
                    if temp not in collection_seq:
                        collection_seq.append(temp)
        return collection_seq


    def _merge_sequences(self, nested_list: list):
        '''
        the prefix of the first sequence is not fixed
        :param nested_list:
        :return:
        '''
        ele_num = len(nested_list)
        if ele_num == 1: return []

        # merge sequence
        merge_seq = nested_list[0]

        for i in range(1, ele_num):
            merge_seq = lists_merge.merge_two_list(merge_seq, nested_list[i])

        # form final sequence
        # remove the elements from the first parent
        final_seq = [[1, nested_list[0][0]]]
        final_seq += merge_seq[1:]
        return [final_seq]

    def generate_sequences(self):
        # get all sequences for each uncovered function
        for ftn_idx in self.ftn_idx_not_covered:
            # get labels and parents
            l_p_dict = self._get_label_parent(ftn_idx)

            if len(l_p_dict) == 0:
                self.ftn_no_sequences.append(ftn_idx)
                continue

            parent_groups = [values for values in l_p_dict.values()]
            all_sequences_ftn = []

            # consider each individual parent
            parents = [p_ftn for group in parent_groups for p_ftn in group]

            if fdg.FDG_global.control_level == 5:
                all_sequences_ftn=self._get_sequences_by_level_5(parents,ftn_idx)
                if len(all_sequences_ftn)>0:
                    self.generated_sequences[ftn_idx] = sorted(all_sequences_ftn, key=len)

                else:self.ftn_no_sequences.append(ftn_idx)
            else:
                # consider each individual parent
                for p_ftn in parents:
                    if p_ftn not in self.all_valid_sequences.keys(): continue
                    p_sequences = self.all_valid_sequences[p_ftn]['sequences']
                    for p_seq in p_sequences:
                        if len(p_seq) == self.fdg.depth_limit:  # consider sequences of length self.fdg.depth_limit
                            if [[self.fdg.depth_limit, p_seq[-1]], ftn_idx] not in all_sequences_ftn:
                                all_sequences_ftn.append([[self.fdg.depth_limit, p_seq[-1]], ftn_idx])

                # consider multiple parents
                if len(parent_groups) > 1:
                    # get all sequences for a function

                    for num_parents in range(1, len(parent_groups)):
                        seq_list=[]
                        if fdg.FDG_global.control_level == 0:  # the prefix of the first sequence is not fixed
                            seq_list=self._get_sequences_by_level_0(parent_groups, num_parents + 1, ftn_idx)

                        elif fdg.FDG_global.control_level <= 2:
                            seq_list= self._get_sequences_by_level_1_2(parent_groups, num_parents + 1, ftn_idx)
                        else:
                            seq_list= self._get_sequences_by_level_3_4(parent_groups, num_parents + 1, ftn_idx)
                        for seq in seq_list:
                            if seq not in all_sequences_ftn:
                                all_sequences_ftn.append(seq)
                # check if there are generated sequences
                if len(all_sequences_ftn) > 0:
                    self.generated_sequences[ftn_idx] = sorted(all_sequences_ftn, key=len)
                else:
                    self.ftn_no_sequences.append(ftn_idx)

        if len(self.generated_sequences.keys()) == 0:
            return False

        # print(f'all sequences generated={self.generated_sequences}')

        # use an 2d array to indicate which sequences are assigned
        ftn_list = []
        ftn_seq_num = []
        for ftn_idx, sequences in self.generated_sequences.items():
            ftn_list.append(ftn_idx)
            ftn_seq_num.append(len(sequences))

        num_ftn = len(ftn_list)
        max_seq_num = max(ftn_seq_num)
        # create a matrix to register all sequences for each uncovered function
        self.generated_sequences_assignment = np.zeros([max_seq_num, num_ftn])
        for idx, ftn_idx in enumerate(ftn_list):
            self.ftn_idx_to_column_idx[ftn_idx] = idx
            self.column_idx_to_ftn_idx[idx] = ftn_idx
            self.generated_sequences_assignment[:, idx][0:ftn_seq_num[idx]] = 1

        return True

    def get_package_dict(self, ftn_not_covered: list):

        # when no sequences are generated
        if len(self.ftn_idx_to_column_idx.keys()) == 0: return {}

        # update: mark all sequences of a function that is covered
        ftn_to_be_assigned = []
        for ftn_idx in self.ftn_idx_to_column_idx.keys():
            if ftn_idx not in ftn_not_covered:
                # unmark all its unassigned sequences
                self.generated_sequences_assignment[:, self.ftn_idx_to_column_idx[ftn_idx]] = 0
            else:
                ftn_to_be_assigned.append(ftn_idx)

        # get sequences for the package
        current_package_dict = {}
        current_package_dict['locate_state'] = []

        # determine package size
        # package_size=len(ftn_to_be_assigned)
        # if self.num_seq_limit>3* len(ftn_to_be_assigned):
        #     package_size=3* len(ftn_to_be_assigned)
        # elif self.num_seq_limit> len(ftn_to_be_assigned):
        #     package_size=self.num_seq_limit
        package_size = 1

        indices = np.where(self.generated_sequences_assignment == 1)
        indices_list = [(x, y) for x, y in zip(indices[0], indices[1])]
        if len(indices_list) > 0:
            # transform sequences into graph
            package_indices = indices_list[0:package_size]
            sequences = []
            for (seq_idx, col_idx) in package_indices:
                # mark assigned sequences 0
                self.generated_sequences_assignment[seq_idx, col_idx] = 0
                # get the sequence
                ftn_idx = self.column_idx_to_ftn_idx[col_idx]
                seq = self.generated_sequences[ftn_idx][seq_idx]
                sequences.append(seq)

            # re-organize sequences in the package
            for seq in sequences:
                if seq[0] not in current_package_dict['locate_state']:
                    current_package_dict['locate_state'] += [seq[0]]
                if 0 not in current_package_dict.keys():
                    current_package_dict[0] = {seq[0][1]: [seq[1]]}
                else:
                    if seq[0][1] not in current_package_dict[0].keys():
                        current_package_dict[0][seq[0][1]] = [seq[1]]
                    else:
                        current_package_dict[0][seq[0][1]] += [seq[1]]

                for i in range(1, len(seq) - 1):
                    if i not in current_package_dict.keys():
                        current_package_dict[i] = {seq[i]: [seq[i + 1]]}
                    else:
                        if seq[i] not in current_package_dict[i].keys():
                            current_package_dict[i][seq[i]] = [seq[i + 1]]
                        else:
                            current_package_dict[i][seq[i]] += [seq[i + 1]]

            # remove repeated elements
            for depth, seq_dict in current_package_dict.items():
                if isinstance(seq_dict, dict):
                    for ftn, values in seq_dict.items():
                        current_package_dict[depth][ftn] = list(set(values))

            self.package_list.append(copy(current_package_dict))
        return current_package_dict

    def get_one_sequence(self, ftn_not_covered: list) -> list:
        # when no sequences are generated
        if len(self.ftn_idx_to_column_idx.keys()) == 0: return []

        # update: mark all sequences of a function that is covered
        for ftn_idx in self.ftn_idx_to_column_idx.keys():
            if ftn_idx not in ftn_not_covered:
                # unmark all its unassigned sequences
                self.generated_sequences_assignment[:, self.ftn_idx_to_column_idx[ftn_idx]] = 0

        package_size = 1
        indices = np.where(self.generated_sequences_assignment == 1)
        indices_list = [(x, y) for x, y in zip(indices[0], indices[1])]
        if len(indices_list) > 0:
            # transform sequences into graph
            package_indices = indices_list[0:package_size]
            sequences = []
            for (seq_idx, col_idx) in package_indices:
                # mark assigned sequences 0
                self.generated_sequences_assignment[seq_idx, col_idx] = 0
                # get the sequence
                ftn_idx = self.column_idx_to_ftn_idx[col_idx]
                seq = self.generated_sequences[ftn_idx][seq_idx]
                return seq

        return []


if __name__ == '__main__':
    # # ftn_info=Function_info('/home/wei/PycharmProjects/Contracts/_wei/Crowdsale.sol', 'Crowdsale')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    # ftn_info = Function_info('/media/sf___share_vms/__contracts_1818/EtherBox.sol', 'EtherBox')
    #
    # function_dict = ftn_info.functions_dict_slither()
    #
    # fdg_object = FDG(function_dict)
    # valid_sequence = {6: [[6], [5, 6]], 2: [[2]], 5: [[5]], 1: [[1], [1, 1], [5, 1]]}
    #
    # seq_object = Sequence(fdg_object, [3,4], valid_sequence, 5)
    # fdg.FDG_global.control_level=2
    # fdg_object.depth_limit=2
    # seq_object.generate_sequences()
    # print(f'sequence={seq_object.generated_sequences}')



    per = permutations([[1, 0], [2, 7], [4, 3]])
    for i in per:
        print(i)

    print(f'===================')
    com = combinations([[1, 0], [2, 7], [4, 3]], r=3)
    for i in com:
        print(i)

    print(f'===================')
    a = [[0, 1], [2, 3], [2, 4]]
    r = [com for com in it.product(*a)]
    print(r)

    print(f'===================')
    print(list(range(2)))

    print(f'===================')
    a = [True, True, False]
    b = np.array([1, 2, 3])
    c = list(np.invert(a))
    print(a)
    print(c)
    print(b[c])

    sequence=[[[2, 6], 3], [[2, 1], 3], [[2, 1], 3], [[1, 6], 1, 3], [[1, 1], 6, 3]]
    unique_seq=[]
    for seq in sequence:
        if seq in unique_seq: continue
        else: unique_seq.append(seq)


    print(f'unique_seq={unique_seq}')

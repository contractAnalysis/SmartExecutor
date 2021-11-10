from copy import copy
from itertools import permutations, combinations
import numpy as np
import itertools as it

import fdg
from fdg import lists_merge, utils

import fdg.FDG_global
"""
sequence generation level:
0:
1: consider parent subsets, shortest sequences, the first sequence is fixed
2: consider parent subsets, shortest sequences, the first sequence is fixed
3:
4:
5: consider all parents once, shortest sequences, the first sequence is not fixed
6: consider all parents once, shortest sequences, topological sorting
7:parent subsets, shortest sequeces, topological sorting
8:parent subsets, shortest sequeces, topological sorting, randomly select m sequences for each target function
9: randomly select parents, shortest sequeces, topological sorting
"""


class Sequence():
    def __init__(self, fdg_=None, ftn_idx_not_covered=[], valid_sequences={}, num_seq_limit=5):
        self.fdg = fdg_
        self.uncovered_ftn_idx = ftn_idx_not_covered  # target function
        self.num_seq_limit = num_seq_limit



        self.valid_sequences_given=valid_sequences
        self.valid_sequences_given_transformed={}  # replace self.all_valid_sequences

        self.assignment_sequence=[]
        self.sequences_generated= []
        self.sequences_generated_cur= {}

        self.uncovered_ftn_idx=ftn_idx_not_covered
        self.uncovered_ftn_idx_waiting=[]

        # transform received valid_sequences:
        for key, value in valid_sequences.items():
            self.valid_sequences_given_transformed[key]=self._valid_sequence_transform(value)





        self.ftn_idx_to_column_idx = {}
        self.column_idx_to_ftn_idx = {}
        self.package_list = []

        self.ftn_no_sequences = []

        
        self.flag_all_ftn_considered=False
        # ==================================
        # prepare data
        if fdg.FDG_global.control_level >=6 and fdg.FDG_global.control_level <=9:
            # remove cycles in edge list
            edges_without_cycles = []
            edges_cycles = []
            for edge in self.fdg.edges_no_f:
                if edge[0] == edge[1]:
                    continue
                if [edge[1], edge[0]] in self.fdg.edges_no_f:
                    edges_cycles.append(edge)
                    continue
                edges_without_cycles.append(edge)

            # get node pairs that violate topological rules
            self.rules = {}
            for edge in edges_without_cycles:
                if edge[1] in self.rules.keys():
                    if edge[0] not in self.rules[edge[1]]:
                        self.rules[edge[1]].append(edge[0])
                else:
                    self.rules[edge[1]] = [edge[0]]

            # save cycles in dict
            self.cycles_dict = {}
            for edge in edges_cycles:
                if edge[0] not in self.cycles_dict.keys():
                    self.cycles_dict[edge[0]] = [edge[1]]
                else:
                    if edge[1] not in self.cycles_dict[edge[0]]:
                        self.cycles_dict[edge[0]].append(edge[1])

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
    def _get_parent(self, ftn_idx) -> dict:
        """
        get parents 
        :param ftn_idx:
        :return:[parent idx]
        """
       
        parents=[]
        sv_read_ary = self.fdg.sv_read[ftn_idx, :]
        sv_read = sv_read_ary[sv_read_ary >= 0]
        for sv in sv_read:
            ftn_write_ary = self.fdg.sv_write[:, sv]
            ftn_write = np.where(ftn_write_ary >= 0)[0]
            parents+=[ftn_i for ftn_i in ftn_write if ftn_i != ftn_idx]
        
        return list(set(parents))
    def _valid_sequence_transform(self,sequence_of_sequences:list):
        if len(sequence_of_sequences) > 0:
            ele_size = [len(item) for item in sequence_of_sequences]
            min_size = min(ele_size)
            shortest_seq = [item for item in sequence_of_sequences if len(item) == min_size]
            return {'sequences': sequence_of_sequences, 'seq_num': len(sequence_of_sequences),
                                                           'shortest': shortest_seq,
                                                           'shortest_depth': min_size}
        else:
            return  {'sequences': [], 'seq_num': 0, 'shortest': [],
                                                           'shortest_depth': 0}

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
                if p_element in self.valid_sequences_given_transformed.keys():
                    p_seq_num.append(self.valid_sequences_given_transformed[p_element]['seq_num'])
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
                    self.valid_sequences_given_transformed[p_ele]['sequences'][index] \
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
                if p_element in self.valid_sequences_given_transformed.keys():
                    p_seq_num.append(len(self.valid_sequences_given_transformed[p_element]['shortest']))
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
                    self.valid_sequences_given_transformed[p_ele]['shortest'][index] \
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
        get one sequence from one parent subset
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

    def _merge_sequences(self, nested_list: list):
        '''
        no permutation, no ordering
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

    # consider parent subsets, prefix of the first sequence is not fixed(no repeated elements in generated sequence)
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
                if p_element in self.valid_sequences_given_transformed.keys():
                    p_seq_num.append(len(self.valid_sequences_given_transformed[p_element]['shortest']))
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
                    self.valid_sequences_given_transformed[p_ele]['shortest'][index] \
                    for p_ele, index in zip(comb, p_seq_index)]
                merge_seq = self._merge_sequences(sequence_list)

                for seq in merge_seq:
                    if len(seq) == 1: continue
                    temp = seq + [ftn_idx]
                    if temp not in collection_seq:
                        collection_seq.append(temp)
        return collection_seq


    # consider all parents at one time, prefix of the first sequence is not fixed
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
        parent_sequences = self._get_parent_sequnces(parents)
        for sequence in parent_sequences:
            permulation_sequence_list= permutations(sequence)
            for seq_list in permulation_sequence_list:
                merge_seq = self._merge_sequences(seq_list)
                for seq in merge_seq:
                    if len(seq) == 1: continue
                    temp = seq + [ftn_idx]
                    if temp not in collection_seq:
                        collection_seq.append(temp)
        return collection_seq

    # topological sorting, consider all parents at one time
    def _get_sequences_by_level_6(self,parents:list, ftn_idx):
        """
        consider all parents
        use shortest sequences to represent parents
        topological sorting in sequence merging process
        :return:
        """
        collection_seq=[]
        parent_sequences = self._get_parent_sequnces(parents)
        for p_seq in parent_sequences:
            sequences = self._get_topological_sequences(p_seq, ftn_idx)
            for seq in sequences:
                if seq not in collection_seq:
                    collection_seq.append(seq)
        return collection_seq

    # topological sorting, consider parent subsets
    def _get_sequences_by_level_7(self,parent_groups: list, num_groups: int, ftn_idx):
        """
        consider parent subsets
        use shortest sequences to represent parents
        topological sorting in sequence merging process
        :return:
        """
        # get sequences through parent combination
        parent_combinations = utils.get_combination(parent_groups, num_groups)

        collection_seq = []  # save the sequences generated
        for comb in parent_combinations:
            parent_sequences = self._get_parent_sequnces(comb)
            for p_seq in parent_sequences:
                sequences = self._get_topological_sequences(p_seq, ftn_idx)
                for seq in sequences:
                    if seq not in collection_seq:
                        collection_seq.append(seq)

        return collection_seq

        # topological sorting, consider parent subsets

    # topological sorting, consider parent subsets, randomly select 5 sequences generated for each function
    def _get_sequences_by_level_8(self, parent_groups: list, num_groups: int, ftn_idx):
        """
        consider parent subsets
        use shortest sequences to represent parents
        topological sorting in sequence merging process
        ramdonly select a specified number of sequences generated for each function
        :return:
        """
        # get sequences through parent combination
        parent_combinations = utils.get_combination(parent_groups, num_groups)

        collection_seq = []  # save the sequences generated
        for comb in parent_combinations:
            parent_sequences = self._get_parent_sequnces(comb)
            for p_seq in parent_sequences:
                sequences=self._get_topological_sequences(p_seq, ftn_idx)
                for seq in sequences:
                    if seq not in collection_seq:
                        collection_seq.append(seq)

        final_sequences = []
        if len(collection_seq)>self.num_seq_limit:
            index=list(range(len(collection_seq)))
            chosen_index = np.random.choice(index, self.num_seq_limit, replace=False)
            final_sequences=[seq for idx, seq in enumerate(collection_seq) if idx in chosen_index]

        return final_sequences

    # topological sorting, randomly select parents multiple times
    def _get_sequences_by_level_9(self, parents:list, ftn_idx):
        """
        randomly select parents, and select multiple times
        use shortest sequences to represent parents
        topological sorting in sequence merging process
        :return:
        """

        # randomly select parents multiple times
        p_len=len(parents)
        p_np_array=np.array(parents)
        p_subsets_selected=[]
        for i in range(self.num_seq_limit):
            selection=[]
            for j in range(p_len):
                selection.append(np.random.choice([True, False]))
            re=list(p_np_array[selection])
            if len(re)==0: # randomly select 2. make sure there is no parent selected
                re=list(np.random.choice(parents,size=2,replace=False))
            p_subsets_selected.append(re)

        collection_seq = []
        # convert each parent subset to sequences
        for p_comb in p_subsets_selected:
            #--------------------------
            # the case only one parent selected
            if len(p_comb)==1:
                collection_seq=self._get_sequence_1_parent_considered(p_comb[0],ftn_idx)
                continue

            # --------------------------
            # the case of more than 1 parent selected
            parent_sequences=self._get_parent_sequnces(p_comb)
            for p_seq in parent_sequences:
                sequences = self._get_topological_sequences(p_seq, ftn_idx)
                for seq in sequences:
                    if seq not in collection_seq:
                        collection_seq.append(seq)
        return collection_seq


    def _get_topological_sequences(self,sequence_list:list,ftn_idx):
        collection_seq=[]
        # get unique nodes
        unique_elements = []
        for seq in sequence_list:
            for s in seq:
                if s not in unique_elements:
                    unique_elements.append(s)

        # combine nodes in a cycle
        super_elements = []
        while len(unique_elements) > 0:
            if unique_elements[0] in self.cycles_dict.keys():
                values_e = self.cycles_dict[unique_elements[0]]
                super_e = []
                for v_e in values_e:
                    if v_e in unique_elements:
                        super_e.append(v_e)
                        unique_elements.remove(v_e)
                super_e.append(unique_elements[0])
                unique_elements.remove(unique_elements[0])

                if super_e not in super_elements:
                    super_elements.append(super_e)
            else:
                super_elements.append([unique_elements[0]])
                unique_elements.remove(unique_elements[0])

        permutation = permutations(super_elements)
        flag = False
        for ele in permutation:
            ele_list = list(ele)
            i = 0
            flag = False
            while i < len(super_elements) - 1:
                rest_elements = [e for k, e_list in enumerate(ele) if k > i for e in e_list]
                not_appear = [self.rules[e] for e in ele_list[i] if e in self.rules.keys()]
                not_appear_list = [e for e_list in not_appear for e in e_list]
                for e in rest_elements:
                    if e in not_appear_list:
                        flag = True
                        break
                if flag:
                    break
                else:
                    i += 1
            if not flag:
                valid_seq = [e for e_list in ele_list for e in e_list]
                valid_seq.append(ftn_idx)
                if valid_seq not in collection_seq:
                    collection_seq.append(valid_seq)

        # form the final sequence, with the first element replaced by a list of length 2:
        # [depth, first element]
        final_sequences = []
        for seq in collection_seq:
            if len(seq) == 1: continue
            first_e = [[1, seq[0]]]
            seq.remove(seq[0])
            final_sequences.append(first_e + seq)
        return final_sequences

    def _get_parent_sequnces(self,parents:list):
        parent_sequences=[]
        p_seq_num = []
        for p_element in parents:
            if p_element in self.valid_sequences_given_transformed.keys():
                p_seq_num.append(len(self.valid_sequences_given_transformed[p_element]['shortest']))
            else:
                p_seq_num.append(0)

        # if one parent does not have sequence, ignore this parent combination
        if p_seq_num.count(0) > 0:
            return []

        # some parent has multiple shortest sequences, so need to do combination
        p_seq_idx_list = [list(range(num)) for num in p_seq_num]
        p_seq_index_comb = [list(com) for com in it.product(*p_seq_idx_list)]
        for p_seq_index in p_seq_index_comb:
            # replace each parent with its shortest sequence
            # (remove 0: constructor is ignored.  reverse: so that parent itself is the last element in its shortest sequence )
            sequence_list = [
                self.valid_sequences_given_transformed[p_ele]['shortest'][index] \
                for p_ele, index in zip(parents, p_seq_index)]
            parent_sequences.append(sequence_list)

        return parent_sequences

    def _get_sequence_1_parent_considered(self,parent_idx:int,ftn_idx:int):
        collection_seq=[]
        if parent_idx not in self.valid_sequences_given_transformed.keys(): return []
        p_sequences = self.valid_sequences_given_transformed[parent_idx]['sequences']
        for p_seq in p_sequences:
            if len(p_seq) == self.fdg.depth_limit:  # consider sequences of length self.fdg.depth_limit
                if [[self.fdg.depth_limit, parent_idx], ftn_idx] not in collection_seq:
                    collection_seq.append([[self.fdg.depth_limit, parent_idx], ftn_idx])
        return collection_seq

    def generate_sequences(self):
        self.sequences_generated_cur= {}
        # delay the sequence generation for functions whose parents are also uncovered.
        if len(self.uncovered_ftn_idx)>0:
            if self.flag_all_ftn_considered:
                return
            if len(self.uncovered_ftn_idx_waiting)>0:
                self.uncovered_ftn_idx=self.uncovered_ftn_idx_waiting
                self.uncovered_ftn_idx_waiting=[]
                
            ftn_to_generate_seq = []            
            for ftn_idx in self.uncovered_ftn_idx:
                parents=self._get_parent(ftn_idx)
                if len(set(parents).intersection(set(self.uncovered_ftn_idx)))>0:
                    self.uncovered_ftn_idx_waiting.append(ftn_idx)
                else:
                    ftn_to_generate_seq.append(ftn_idx)
                    
            if len(ftn_to_generate_seq)==0:# 
                self.uncovered_ftn_idx=self.uncovered_ftn_idx_waiting
                self.uncovered_ftn_idx_waiting=[]
                self.flag_all_ftn_considered=True
            else:
                self.uncovered_ftn_idx=ftn_to_generate_seq
                if len(self.uncovered_ftn_idx_waiting)==0:
                    self.flag_all_ftn_considered=True
                
        else: return # no need to generate sequence
            
        
        
        #==================================
        # get all sequences for each uncovered function
        for ftn_idx in self.uncovered_ftn_idx:
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
                # consider all parents(shortest sequences to represent parents)
                # the first sequence is not fixed
                all_sequences_ftn=self._get_sequences_by_level_5(parents,ftn_idx)
                if len(all_sequences_ftn)>0:
                    self.sequences_generated_cur[ftn_idx] = sorted(all_sequences_ftn, key=len)
                else:self.ftn_no_sequences.append(ftn_idx)
                
            elif fdg.FDG_global.control_level==6:
                # consider all parents(shortest sequences to represent parents)
                # consider topological sorting in merged sequences
                all_sequences_ftn = self._get_sequences_by_level_6(parents, ftn_idx)
                if len(all_sequences_ftn) > 0:
                    self.sequences_generated_cur[ftn_idx] = sorted(all_sequences_ftn, key=len)

                else:
                    self.ftn_no_sequences.append(ftn_idx)
            elif fdg.FDG_global.control_level==9:
                # consider all parents(shortest sequences to represent parents)
                # consider topological sorting in merged sequences
                all_sequences_ftn = self._get_sequences_by_level_9(parents, ftn_idx)
                if len(all_sequences_ftn) > 0:
                    self.sequences_generated_cur[ftn_idx] = sorted(all_sequences_ftn, key=len)

                else:
                    self.ftn_no_sequences.append(ftn_idx)
            else:
                #========================================
                # consider parent subsets
                # case 1: consider each individual parent
                for p_ftn in parents:
                    if p_ftn not in self.valid_sequences_given_transformed.keys(): continue
                    p_sequences = self.valid_sequences_given_transformed[p_ftn]['sequences']
                    for p_seq in p_sequences:
                        if len(p_seq) == self.fdg.depth_limit:  # consider sequences of length self.fdg.depth_limit
                            if [[self.fdg.depth_limit, p_seq[-1]], ftn_idx] not in all_sequences_ftn:
                                all_sequences_ftn.append([[self.fdg.depth_limit, p_seq[-1]], ftn_idx])

                # case 2: consider multiple parents
                if len(parent_groups) > 1:
                    # get all sequences for a function
                    for num_parents in range(1, len(parent_groups)):
                        seq_list=[]
                        if fdg.FDG_global.control_level == 0:  # the prefix of the first sequence is not fixed
                            seq_list=self._get_sequences_by_level_0(parent_groups, num_parents + 1, ftn_idx)

                        elif fdg.FDG_global.control_level <= 2:
                            seq_list= self._get_sequences_by_level_1_2(parent_groups, num_parents + 1, ftn_idx)
                        elif fdg.FDG_global.control_level<=4:
                            seq_list= self._get_sequences_by_level_3_4(parent_groups, num_parents + 1, ftn_idx)
                        elif fdg.FDG_global.control_level==7:
                            seq_list=self._get_sequences_by_level_7(parent_groups,num_parents+1,ftn_idx)
                        elif fdg.FDG_global.control_level == 8:
                            seq_list = self._get_sequences_by_level_8(parent_groups, num_parents + 1, ftn_idx)

                        for seq in seq_list:
                            if seq not in all_sequences_ftn:
                                all_sequences_ftn.append(seq)
                # check if there are generated sequences
                if len(all_sequences_ftn) > 0:
                    self.sequences_generated_cur[ftn_idx] = sorted(all_sequences_ftn, key=len)
                else:
                    self.ftn_no_sequences.append(ftn_idx)

        if len(self.sequences_generated_cur.keys()) == 0:
            return

        # print(f'all sequences generated={self.generated_sequences}')
        self.sequences_generated.append(self.sequences_generated_cur)

        # use an 2d array to indicate which sequences are assigned
        ftn_list = []
        ftn_seq_num = []
        for ftn_idx, sequences in self.sequences_generated_cur.items():
            ftn_list.append(ftn_idx)
            ftn_seq_num.append(len(sequences))

        num_ftn = len(ftn_list)
        max_seq_num = max(ftn_seq_num)
        # create a matrix to register all sequences for each uncovered function
        self.assignment_sequence= np.zeros([max_seq_num, num_ftn])
        for idx, ftn_idx in enumerate(ftn_list):
            self.ftn_idx_to_column_idx[ftn_idx] = idx
            self.column_idx_to_ftn_idx[idx] = ftn_idx
            self.assignment_sequence[:, idx][0:ftn_seq_num[idx]] = 1
        return



    def get_one_sequence(self, ftn_not_covered: list) -> list:
        sequences=[]
        if len(self.assignment_sequence)>0:
            # update self.assignment_sequence
            for ftn_idx in self.uncovered_ftn_idx:
                if ftn_idx not in ftn_not_covered:
                    # unmark all its unassigned sequences
                    self.assignment_sequence[:, self.ftn_idx_to_column_idx[ftn_idx]] = 0

            # there are still generated sequences not assigned yet
            if 1 in self.assignment_sequence:
                # get sequences, the package size is 1
                sequences=self.assign_mark_a_sequence(1)[0]
                return sequences


        # update valid_sequences_given_transformed
        for key in self.valid_sequences_given.keys():
            if key not in self.valid_sequences_given_transformed.keys():
                value=self.valid_sequences_given[key]
                self.valid_sequences_given_transformed[key]=self._valid_sequence_transform(value)

        # generate and assigne sequences
        self.generate_sequences()
        if len(self.sequences_generated_cur) == 0:
            return []
        sequences = self.assign_mark_a_sequence(1)[0]
        return sequences

    def  assign_mark_a_sequence(self,package_size:int):

        indices = np.where(self.assignment_sequence== 1)
        indices_list = [(x, y) for x, y in zip(indices[0], indices[1])]
        sequences = []
        if len(indices_list) > 0:
            # get a list of sequences
            package_indices = indices_list[0:package_size]
            for (seq_idx, col_idx) in package_indices:
                # mark assigned sequences 0
                self.assignment_sequence[seq_idx, col_idx] = 0
                # get the sequence
                ftn_idx = self.column_idx_to_ftn_idx[col_idx]
                seq = self.sequences_generated_cur[ftn_idx][seq_idx]
                sequences.append(seq)
        return sequences




    # def get_one_sequence(self, ftn_not_covered: list) -> list:
    #     # when no sequences are generated
    #     if len(self.ftn_idx_to_column_idx.keys()) == 0: return []
    #
    #     # update: mark all sequences of a function that is covered
    #     for ftn_idx in self.ftn_idx_to_column_idx.keys():
    #         if ftn_idx not in ftn_not_covered:
    #             # unmark all its unassigned sequences
    #             self.generated_sequences_assignment[:, self.ftn_idx_to_column_idx[ftn_idx]] = 0
    #
    #     package_size = 1
    #     indices = np.where(self.generated_sequences_assignment == 1)
    #     indices_list = [(x, y) for x, y in zip(indices[0], indices[1])]
    #     if len(indices_list) > 0:
    #         # transform sequences into graph
    #         package_indices = indices_list[0:package_size]
    #         sequences = []
    #         for (seq_idx, col_idx) in package_indices:
    #             # mark assigned sequences 0
    #             self.generated_sequences_assignment[seq_idx, col_idx] = 0
    #             # get the sequence
    #             ftn_idx = self.column_idx_to_ftn_idx[col_idx]
    #             seq = self.generated_sequences[ftn_idx][seq_idx]
    #             return seq
    #
    #     return []




def _merge_sequences(nested_list:list):
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

    # # form final sequence
    # # remove the elements from the first parent
    # final_seq = [[1, nested_list[0][0]]]
    # final_seq += merge_seq[1:]
    return [merge_seq]


if __name__ == '__main__':
    # # ftn_info=Function_info('/home/wei/PycharmProjects/Contracts/_wei/Crowdsale.sol', 'Crowdsale')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    # # ftn_info = Function_info('/media/sf___share_vms/__contracts_1818/EtherBox.sol', 'EtherBox')
    #
    # function_dict = ftn_info.functions_dict_slither()
    #
    # fdg_object = FDG(function_dict)
    # # valid_sequence = {6: [[6], [5, 6]], 2: [[2]], 5: [[5]], 1: [[1], [1, 1], [5, 1]]}
    # valid_sequence={13: {'sequences': [[13], [2, 13], [2, 2, 13]], 'seq_num': 3, 'shortest': [[13]], 'shortest_depth': 1}, 4: {'sequences': [[4]], 'seq_num': 1, 'shortest': [[4]], 'shortest_depth': 1}, 10: {'sequences': [[10], [2, 10], [2, 2, 10]], 'seq_num': 3, 'shortest': [[10]], 'shortest_depth': 1}, 7: {'sequences': [[7]], 'seq_num': 1, 'shortest': [[7]], 'shortest_depth': 1}, 2: {'sequences': [[2], [2, 2], [2, 2, 2]], 'seq_num': 3, 'shortest': [[2]], 'shortest_depth': 1}, 9: {'sequences': [[9]], 'seq_num': 1, 'shortest': [[9]], 'shortest_depth': 1}, 8: {'sequences': [[8]], 'seq_num': 1, 'shortest': [[8]], 'shortest_depth': 1}, 12: {'sequences': [[10, 12], [2, 10, 12]], 'seq_num': 2, 'shortest': [[10, 12]], 'shortest_depth': 2}, 11: {'sequences': [[10, 11], [2, 10, 11], [10, 11, 11]], 'seq_num': 3, 'shortest': [[10, 11]], 'shortest_depth': 2}, 6: {'sequences': [[10, 12, 6]], 'seq_num': 1, 'shortest': [[10, 12, 6]], 'shortest_depth': 3}, 3: {'sequences': [[10, 12, 3]], 'seq_num': 1, 'shortest': [[10, 12, 3]], 'shortest_depth': 3}, 5: {'sequences': [[10, 12, 5]], 'seq_num': 1, 'shortest': [[10, 12, 5]], 'shortest_depth': 3}}
    # seq_object = Sequence(fdg_object, [3,4], valid_sequence, 5)
    # fdg.FDG_global.control_level=2
    # fdg_object.depth_limit=2
    # seq_object.generate_sequences()
    # print(f'sequence={seq_object.generated_sequences}')
    print(np.random.choice([1,2,3,4],size=2))
    selection=[]
    targets=[2,4,5,7,8]
    for i in range(6):
        selection=[]
        for i in range(5):
            selection.append(np.random.choice([True,False]))

        print(f'randomly choose: {selection}')
        print(f'randomly choose: {list(np.array(targets)[selection])}')

    index = list(range(10))
    chosen_index = np.random.choice(index,5,replace=False)
    print(index)
    print(chosen_index)

    #
    # edges = [[2, 2], [2, 10], [2, 13], [3, 3], [3, 5], [3, 11], [3, 14], [5, 5], [5, 3], [5, 11], [5, 14], [6, 5],
    #          [8, 5], [9, 5], [10, 11], [10, 12], [11, 3], [11, 5], [11, 11], [11, 14], [12, 3], [12, 5], [12, 6],
    #          [12, 11], [13, 14], [14, 3], [14, 5], [14, 11], [14, 14]]
    #
    # rules={}
    #
    # sequences_list = [[10,12,3], [10,12,5], [10,11], [13]]
    # unique_elements=[]
    #
    # # remove cycles in edge list
    # edges_without_cycles=[]
    # edges_cycles=[]
    # for edge in edges:
    #     if edge[0]==edge[1]:
    #         continue
    #     if [edge[1],edge[0]] in edges:
    #         edges_cycles.append(edge)
    #         continue
    #     edges_without_cycles.append(edge)
    #
    # # get node pairs that violate topological rules
    # for edge in edges_without_cycles:
    #     if edge[1] in rules.keys():
    #         if edge[0] not in rules[edge[1]]:
    #             rules[edge[1]].append(edge[0])
    #     else: rules[edge[1]]=[edge[0]]
    #
    # # save cycles in dict
    # cycles_dict={}
    # for edge in edges_cycles:
    #     if edge[0] not in cycles_dict.keys():
    #         cycles_dict[edge[0]]=[edge[1]]
    #     else:
    #         if edge[1] not in cycles_dict[edge[0]]:
    #             cycles_dict[edge[0]].append(edge[1])
    #
    #
    # # get unique nodes
    # for seq in sequences_list:
    #     for s in seq:
    #         if s not in unique_elements:
    #             unique_elements.append(s)
    #
    # # combine nodes in a cycle
    # super_elements=[]
    # e_temp=unique_elements
    # while len(unique_elements)>0:
    #     if unique_elements[0] in cycles_dict.keys():
    #         values_e=cycles_dict[unique_elements[0]]
    #         super_e=[]
    #         for v_e in values_e:
    #             if v_e in unique_elements:
    #                 super_e.append(v_e)
    #                 unique_elements.remove(v_e)
    #         super_e.append(unique_elements[0])
    #         unique_elements.remove(unique_elements[0])
    #
    #         if super_e not in super_elements:
    #             super_elements.append(super_e)
    #     else:
    #         super_elements.append([unique_elements[0]])
    #         unique_elements.remove(unique_elements[0])
    #
    #
    # sequences=[]
    # permutation=permutations(super_elements)
    # flag=False
    # for ele in permutation:
    #     ele_list=list(ele)
    #     i=0
    #     flag=False
    #     while i<len(super_elements)-1:
    #         rest_elements=[e for k, e_list in enumerate(ele) if k>i for e in e_list]
    #         not_appear=[rules[e] for e in ele_list[i] if e in rules.keys() ]
    #         not_appear_list=[e for e_list in not_appear for e in e_list]
    #         for e in rest_elements:
    #             if e in not_appear_list:
    #                 flag=True
    #                 break
    #         if flag:
    #             break
    #         else:
    #             i+=1
    #     if not flag:
    #         sequences.append([e for e_list in ele_list for e in e_list])
    #
    #
    # i=0
    # for ele in sequences:
    #     print(f'{i}: {ele}')
    #     i+=1



    # per = permutations([[1, 0], [2, 7], [4, 3]])
    # for i in per:
    #     print(i)
    #
    # print(f'===================')
    # com = combinations([[1, 0], [2, 7], [4, 3]], r=3)
    # for i in com:
    #     print(i)
    #
    # print(f'===================')
    # a = [[0, 1], [2, 3], [2, 4]]
    # r = [com for com in it.product(*a)]
    # print(r)
    #
    # print(f'===================')
    # print(list(range(2)))
    #
    # print(f'===================')
    # a = [True, True, False]
    # b = np.array([1, 2, 3])
    # c = list(np.invert(a))
    # print(a)
    # print(c)
    # print(b[c])
    #
    # sequence=[[[2, 6], 3], [[2, 1], 3], [[2, 1], 3], [[1, 6], 1, 3], [[1, 1], 6, 3]]
    # unique_seq=[]
    # for seq in sequence:
    #     if seq in unique_seq: continue
    #     else: unique_seq.append(seq)
    # print(f'unique_seq={unique_seq}')

    holoToken_generated_sequences=[[[1, 10], 12, 3, 5, 11, 13, 14], [[1, 10], 12, 3, 5, 13, 11, 14], [[1, 10], 12, 3, 11, 5, 13, 14], [[1, 10], 12, 3, 11, 13, 5, 14], [[1, 10], 12, 3, 13, 5, 11, 14], [[1, 10], 12, 3, 13, 11, 5, 14], [[1, 10], 12, 5, 3, 11, 13, 14], [[1, 10], 12, 5, 3, 13, 11, 14], [[1, 10], 12, 5, 11, 3, 13, 14], [[1, 10], 12, 5, 11, 13, 3, 14], [[1, 10], 12, 5, 13, 3, 11, 14], [[1, 10], 12, 5, 13, 11, 3, 14], [[1, 10], 11, 12, 3, 5, 13, 14], [[1, 10], 11, 12, 3, 13, 5, 14], [[1, 10], 11, 12, 5, 3, 13, 14], [[1, 10], 11, 12, 5, 13, 3, 14], [[1, 10], 11, 13, 12, 3, 5, 14], [[1, 10], 11, 13, 12, 5, 3, 14], [[1, 13], 10, 12, 3, 5, 11, 14], [[1, 13], 10, 12, 3, 11, 5, 14], [[1, 13], 10, 12, 5, 3, 11, 14], [[1, 13], 10, 12, 5, 11, 3, 14], [[1, 13], 10, 11, 12, 3, 5, 14], [[1, 13], 10, 11, 12, 5, 3, 14]]

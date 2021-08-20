from copy import copy
from itertools import permutations,combinations
import numpy as np
import itertools as it
from fdg import lists_merge



class Sequence():
    def __init__(self, fdg=None,ftn_idx_not_covered=[],num_seq_limit=5):
        self.fdg=fdg
        self.ftn_idx_not_covered=ftn_idx_not_covered# target function
        self.num_seq_limit=num_seq_limit
        self.ftn_seq_and_shortest_seq_dict={}# for save all sequences for each function from updated FDG


        self.all_sequences = {}
        self.all_sequences_assignment=[]
        self.ftn_list_assignment=[]
        self.package_list=[]

        self.shortest_flag=True


    def _get_label_parent(self,ftn_idx)->dict:
        """
        get parents and their labels on the orginal fdg
        :param ftn_idx:
        :return:dict: key: label index, value: parent indices
        """
        l_p_dict={}
        sv_read_ary=self.fdg.sv_read[ftn_idx,:]
        sv_read=sv_read_ary[sv_read_ary>=0]
        for sv in sv_read:
            ftn_write_ary=self.fdg.sv_write[:,sv]
            ftn_write=np.where(ftn_write_ary>=0)[0]
            l_p_dict[sv]=[ftn_i for ftn_i in ftn_write if ftn_i!=ftn_idx]
        return l_p_dict

    def _generate_backward_graph(self,new_fdg:bool)->dict:
        """
        build graph backward
        from lead nodes to contructor
        :return:
        """
        graph_bk={}
        if new_fdg:
            fdg_3d_array=self.fdg.fdg_3d_array_new
        else:fdg_3d_array=self.fdg.fdg_3d_array

        # graph_backward_parent_label={}
        for i in range(fdg_3d_array.shape[0]-1,-1,-1):
            ftn_reached=self.fdg.ftn_reached_at_a_depth(i,new_fdg)
            for ftn_to in ftn_reached:
                ftn_parents_all=fdg_3d_array[i,ftn_to,:]
                ftn_parents=list(np.where(ftn_parents_all>=0)[0])
                if len(ftn_parents)==0:continue
                if ftn_to in ftn_parents:
                    ftn_parents.remove(ftn_to)
                if ftn_to not in graph_bk.keys():
                    graph_bk[ftn_to]=ftn_parents
                else:
                    graph_bk[ftn_to]+=ftn_parents
        # remove repeated elements
        for key in graph_bk.keys():
            graph_bk[key]=list(set(graph_bk[key]))
        if new_fdg:
            self.graph_new=graph_bk
        else:self.graph_original=graph_bk

    def _generate_sequences(self,ftn_idx:int)->dict:
        """
        generate sequences from both original and new fdg
        :param ftn_idx:
        :return:{'original':[...],'new':[...]}
        """
        assert self.graph_original
        assert self.graph_new

        def _get_sequences(graph, function_idx, visited, path, seq):
            """
            The general method to get seqences from function_idx to leaf nodes in the graph
            :param graph: a dict
            :param function_idx: the source
            :param visited: [# of nodes]
            :param path:[]
            :param seq: save the resulted sequences

            """
            # checking all the visited nodes
            visited[function_idx] = True
            path.append(function_idx)

            if function_idx not in graph:  # leaf nodes
                # se.append(list(np.copy(path))[1::][::-1])
                seq.append(list(np.copy(path)))
            else:
                for i in graph[function_idx]:
                    if visited[i] == False: 
                        _get_sequences(graph, i, visited, path, seq)
            path.pop()
            visited[function_idx] = False

        # get sequences from original FDG
        seq_ori=[]
        # path=[]
        # visited = [False] * self.fdg.num_ftn
        # _get_sequences(self.graph_original,ftn_idx,visited,path,seq_ori)

        seq_new = []
        # get sequences from new FDG (some functions may not appear in FDG)
        if ftn_idx in self.graph_new:
            seq_new = []
            path=[]
            visited = [False] * self.fdg.num_ftn
            _get_sequences(self.graph_new, ftn_idx, visited, path, seq_new)

        return {'original':seq_ori,'new':seq_new}


    def _get_sequences_and_shortest_sequences(self):
        """
        get sequences and shortest sequence from new FDG
        ignore sequences from original FDG
        :return:
        """
        ftn_seq_dict = {}
        for ftn_idx in range(1, self.fdg.num_ftn):
            seqs = {}
            # get sequences
            sequences = self._generate_sequences(ftn_idx)

            # only consider sequences from new FDG
            for key, value in sequences.items():
                if key=='new':
                    if len(value) > 0:
                        ele_size = [len(item) for item in value]
                        min_size = min(ele_size)
                        shortest_seq = [item for item in value if len(item) == min_size]
                        seqs[key] = {'sequences': value, 'seq_num':len(value),'shortest':shortest_seq,'shortest_depth':min_size-1}

                    else:
                        seqs[key] = {'sequences': [], 'seq_num':0, 'shortest':[],'shortest_depth':0}

            ftn_seq_dict[ftn_idx] = seqs
        return ftn_seq_dict

    def _get_combination(self, list_for_comb):
            """
            -group indices with the same value [[1,4], [2,6],[5]]
            -from each group, choose 1 index
            -generate combinations with length from 2 to the number of groups
            :param list_for_comb: like [-1, 0, 1, -1, 0, 2, 1]
            :return: [(1, 2), (1, 6), (4, 2), (4, 6), (1, 5), (4, 5), (2, 5), (6, 5), (1, 2, 5), (1, 6, 5), (4, 2, 5), (4, 6, 5)]

            """

            com_re = []
            # do combination from length 2 to the number of groups
            num_groups = len(list_for_comb)
            group_indices = range(num_groups)  # used for groups selection for combination

            for c_len in range(2, num_groups + 1, 1):

                # select c_len groups for combinations
                com_groups = it.combinations(group_indices, c_len)
                # for each c_len-group combination, get combination
                for idx, g_indices in enumerate(com_groups):
                    # combine groups to form a list, with each elment is the group
                    g_for_com = [list_for_comb[g_indices[i]] for i in range(len(g_indices))]
                    # do combination
                    com_re += [com for com in it.product(*g_for_com)]

            return com_re

    def _get_combination(self, list_for_comb,comb_length:int):
        """

        :param list_for_comb: [[1,4], [2,6],[5]]
        :param comb_length: 2 (two elements in a combination)
        :return: [(1, 2), (1, 6), (4, 2), (4, 6), (1, 5), (4, 5), (2, 5), (6, 5)]
        """
        com_re = []
        # do combination with length 2
        num_groups = len(list_for_comb)
        if num_groups<comb_length:return []

        group_indices = range(num_groups)  # used for groups selection for combination

        # select comb_length groups for combinations
        com_groups = it.combinations(group_indices, comb_length)

        # for each group combination, get combination
        for idx, g_indices in enumerate(com_groups):
            # combine groups to form a list, with each elment is the group
            g_for_com = [list_for_comb[g_indices[i]] for i in range(len(g_indices))]
            # do combination
            com_re += [com for com in it.product(*g_for_com)]

        return com_re



    """
    updated sequence organization
    generate and execute sequences level by level
    level: reflect the number of parents considerd to generate sequences for a target function
    """
    def prepare_sequence_generation(self):
        # get shortest sequences for all functions except constructor
        self._generate_backward_graph(True)  # build graph based on new fdg
        self._generate_backward_graph(False)  # build graph based on origial fdg

        # key-value format: target function: [[target function, other function, constructor],...]
        self.ftn_seq_and_shortest_seq_dict = self._get_sequences_and_shortest_sequences()

    # def _get_sequences_by_level(self, parent_groups:list,level:int,ftn_idx):
    #     """
    #   consider all sequences
    #     :param parent_groups:[[2,3,5],[3,6]]
    #     :param level : number of parents to consider
    #     :param ftn_idx : the child
    #     :return:
    #     """
    #     assert self.ftn_seq_and_shortest_seq_dict
    #
    #     # get sequences through parent combination
    #     parent_combinations = self._get_combination(parent_groups, level)
    #
    #     collection_seq = []  # save the sequences generated
    #     for comb in parent_combinations:
    #         # check if each parent has shortest sequence or not
    #         p_seq_num = [self.ftn_seq_and_shortest_seq_dict[p_element]['new']['seq_num'] for p_element in
    #                      comb]
    #         # if one parent does not have sequence, ignore this parent combination
    #         if p_seq_num.count(0) > 0: continue
    #
    #         # some parent has multiple shortest sequences, so need to do combination
    #         p_seq_num_list = [list(range(num)) for num in p_seq_num]
    #         p_seq_index_comb = [list(com) for com in it.product(*p_seq_num_list)]
    #         for p_seq_index in p_seq_index_comb:
    #             # replace each parent with its shortest sequence
    #             # (remove 0: constructor is ignored.  reverse: so that parent itself is the last element in its shortest sequence )
    #             nested_sequence = [
    #                 self.ftn_seq_and_shortest_seq_dict[p_ele]['new']['sequences'][index][0:-1][::-1] \
    #                 for p_ele, index in zip(comb, p_seq_index)]
    #             merge_seq = self._merge_sequences(nested_sequence)
    #

    #             for seq in merge_seq:
    #                 seq.append(ftn_idx)
    #                 if seq not in collection_seq:
    #                     collection_seq.append(merge_seq)
    #
    #     return collection_seq



    def _get_sequences_by_level(self, parent_groups:list,level:int,ftn_idx):
        """
        Only consider shortest sequences
        :param parent_groups:[[2,3,5],[3,6]]
        :param level : number of parents to consider
        :param ftn_idx : the child
        :return:
        """
        assert self.ftn_seq_and_shortest_seq_dict

        # get sequences through parent combination
        parent_combinations = self._get_combination(parent_groups, level)

        collection_seq = []  # save the sequences generated
        for comb in parent_combinations:
            # check if each parent has shortest sequence or not
            p_seq_num = [len(self.ftn_seq_and_shortest_seq_dict[p_element]['new']['shortest'] )for p_element in
                         comb]
            # if one parent does not have sequence, ignore this parent combination
            if p_seq_num.count(0) > 0: continue

            # some parent has multiple shortest sequences, so need to do combination
            p_seq_num_list = [list(range(num)) for num in p_seq_num]
            p_seq_index_comb = [list(com) for com in it.product(*p_seq_num_list)]
            for p_seq_index in p_seq_index_comb:
                # replace each parent with its shortest sequence
                # (remove 0: constructor is ignored.  reverse: so that parent itself is the last element in its shortest sequence )
                nested_sequence = [
                    self.ftn_seq_and_shortest_seq_dict[p_ele]['new']['shortest'][index][0:-1][::-1] \
                    for p_ele, index in zip(comb, p_seq_index)]

                merge_seq = self._merge_sequences_ordered(nested_sequence)
                #merge_seq = self._merge_sequences_all_permutation(nested_sequence)

                for seq in merge_seq:
                    temp=seq+[ftn_idx]
                    if temp not in collection_seq:
                        collection_seq.append(temp)

        return collection_seq

    def _merge_sequences_all_permutation(self, nested_list: list):
        """
        each permulation of nested list is created one sequence
        :param nested_list:
        :return:
        """
        result=[]
        ele_num = len(nested_list)
        if ele_num == 1: return []

        permulation_nested_list=permutations(nested_list)
        for per_tuple in permulation_nested_list:
            temp_list=list(per_tuple)
            p_first = temp_list[0]
            p_first_length = len(temp_list[0])
            merge_seq = p_first
            for i in range(1,ele_num):
                merge_seq = lists_merge.merge_fix_list_1_specified_lenth(merge_seq, temp_list[i],
                                                                         p_first_length)
            # form final sequence
            # remove the elements from the first parent
            final_seq = [[p_first_length, p_first[-1]]]
            final_seq += merge_seq[p_first_length:]
            if final_seq not in result:
                result.append(final_seq)
        return result

    def _merge_sequences_ordered(self,nested_list:list):
            ele_num=len(nested_list)
            if ele_num==1: return []
            # order elements in nested_list based on len
            order_nested_list=sorted(nested_list,key=len,reverse=True)

            # merge sequence
            p_first = order_nested_list[0]
            p_first_length = len(order_nested_list[0])
            merge_seq = p_first

            for i in range(1, ele_num):
                merge_seq = lists_merge.merge_fix_list_1_specified_lenth(merge_seq, order_nested_list[i],
                                                                         p_first_length)

            # form final sequence
            # remove the elements from the first parent
            final_seq=[[p_first_length, p_first[-1]]]
            final_seq+=merge_seq[p_first_length:]
            return  [final_seq]


    def get_all_sequences(self):
        assert self.ftn_seq_and_shortest_seq_dict
        # get all sequences for each function
        for ftn_idx in self.ftn_idx_not_covered:
            # get labels and parents
            l_p_dict = self._get_label_parent(ftn_idx)
            if len(l_p_dict) == 0: continue

            parent_groups = [values for values in l_p_dict.values()]
            if len(parent_groups) <=1: continue

            # get all sequences for a function
            all_sequences_ftn=[]
            for num_parents in range(1,len(parent_groups)):
                all_sequences_ftn+=self._get_sequences_by_level(parent_groups,num_parents+1,ftn_idx)

            if len(all_sequences_ftn)>0:
                self.all_sequences[ftn_idx]=sorted(all_sequences_ftn,key=len)

        if len(self.all_sequences.keys())==0:return False, {}

        # use an 2d array to indicate which sequences are assigned
        ftn_list = []
        ftn_seq_num = []
        for ftn_idx, sequences in self.all_sequences.items():
            ftn_list.append(ftn_idx)
            ftn_seq_num.append(len(sequences))

        self.ftn_list_assignment = ftn_list
        num_ftn = len(ftn_list)
        max_seq_num = max(ftn_seq_num)
        self.all_sequences_assignment= np.zeros([max_seq_num, num_ftn])
        for idx in range(num_ftn):
            self.all_sequences_assignment[:, idx][0:ftn_seq_num[idx]] = 1


        # return functions that do not have sequences
        return True, list(set(self.ftn_idx_not_covered).difference(set(self.all_sequences.keys())))


    def get_package_dict(self,ftn_not_covered:list):
        if len(self.ftn_list_assignment)==0:return {}

        # update: mark all sequences of a function that is covered
        num_ftn_involved=[]
        for i,ftn in enumerate(self.ftn_list_assignment):
            if ftn not in ftn_not_covered:
                self.all_sequences_assignment[:,i]=0
            else:num_ftn_involved.append(ftn)

        # get sequences for the package
        current_package_dict={}
        current_package_dict['locate_state'] = []

        package_size=self.num_seq_limit if self.num_seq_limit>len(num_ftn_involved) else len(num_ftn_involved)
        indices=np.where(self.all_sequences_assignment==1)
        indices_list=[(x,y) for x,y in zip(indices[0],indices[1])]
        if len(indices_list)>0:
            # transform sequences into graph
            package_indices=indices_list[0:package_size]
            sequences=[]
            for (seq_idx,col_idx) in package_indices:
                ftn_idx=self.ftn_list_assignment[col_idx]
                seq=self.all_sequences[ftn_idx][seq_idx]
                sequences.append(seq)

            for seq in sequences:
                if seq[0] not in current_package_dict['locate_state']:
                    current_package_dict['locate_state'] += [seq[0]]
                if 0 not in current_package_dict.keys():
                    current_package_dict[0] = {seq[0][1]: [seq[1]]}
                else:
                    if seq[0][1] not in current_package_dict[0].keys():
                        current_package_dict[0][seq[0][1]]=[seq[1]]
                    else:
                        current_package_dict[0][seq[0][1]]+=[seq[1]]

                for i in range(1,len(seq)-1):
                    if i not in current_package_dict.keys():
                        current_package_dict[i]={seq[i]:[seq[i+1]]}
                    else:
                        if seq[i] not in current_package_dict[i].keys():
                            current_package_dict[i][seq[i]] = [seq[1+1]]
                        else:
                            current_package_dict[i][seq[i]] += [seq[1+1]]

            # mark assigned sequences 0
            for row_idx,col_idx in package_indices:
                self.all_sequences_assignment[row_idx,col_idx]=0

            # remove repeated elements
            for depth, seq_dict in current_package_dict.items():
                if isinstance(seq_dict,dict):
                    for ftn,values in seq_dict.items():
                        current_package_dict[depth][ftn]=list(set(values))

            self.package_list.append(copy(current_package_dict))
        return current_package_dict














if __name__=='__main__':
    # # ftn_info=Function_info('/home/wei/PycharmProjects/Contracts/_wei/Crowdsale.sol', 'Crowdsale')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    #
    # function_dict=ftn_info.functions_dict_slither()
    #
    # fdg = FDG(function_dict)
    # function_list = [1, 5, 7, 6, 10]
    # fdg.build_fdg_3d_array(function_list)

    per=permutations([[1,0],[2,7],[4,3]])
    for i in per:
        print(i)

    print(f'===================')
    com=combinations([[1,0],[2,7],[4,3]],r=1)
    for i in com:
        print(i)
    # r={11: [[(3, 2), 10, 11], [(3, 3), 10, 11], [(2, 8), 10, 11]], 3: [[(3, 2), 5, 3], [(2, 8), 5, 3], [(3, 2), 6, 3], [(2, 8), 6, 3], [(2, 9), 5, 3]]}
    # for k,v in r.items():
    #
    #         r[k]+=[20]
    #
    # print(r)
    # r={1:{2:[2],4:[7]}}
    # for i in range(15):
    #     if i in r[1].keys():
    #         print(f'{i} is in r.')
    # a=[0,0,0]
    # if a== [0] * 3:
    #     print('true')

    # path = []
    # visited= [False] * fdg.num_ftn
    # seq_s = []
    # seq = Sequence(fdg, 'withdraw()')
    # seq.get_sequences(fdg.graph_dict_readwrite, 3, visited, path, seq_s)
    # print(seq_s)



    # seq=[[1],[2,3,1],[2],[1,2],[5]]
    # seq_size = [len(item) for item in seq]
    # min_size = min(seq_size)
    # a=[item for i, item in enumerate(seq) if len(item)==min_size]
    # print(a)

    # seq= Sequence()
    # re=seq._get_combination([[1,4], [2,6],[5]],2)
    # print(re)


    a=[[0,1],[2,3],[2,4]]
    r=[com for com in it.product(*a)]
    print(r)
    print(list(range(2)))

    pass









import numpy as np
import itertools as it

from fdg.FDG_3d_array import FDG
from fdg.funtion_info import Function_info


class Sequence():
    def __init__(self, fdg=None,ftn_idx_not_covered=[],num_seq_limit=5):
        self.fdg=fdg
        self.ftn_idx_not_covered=ftn_idx_not_covered# target function
        self.num_seq_limit=num_seq_limit
        self.ftn_seq_and_shortest_seq_dict={}

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
        graph_backward={}
        if new_fdg:
            fdg_3d_array=self.fdg.fdg_3d_array_new
        else:
            fdg_3d_array=self.fdg.fdg_3d_array

        # graph_backward_parent_label={}
        for i in range(fdg_3d_array.shape[0]-1,-1,-1):
            ftn_reached=self.fdg.ftn_reached_at_a_depth(i,new_fdg)
            for ftn_from in ftn_reached:
                ftn_parents_all=fdg_3d_array[i,ftn_from,:]
                ftn_parents=np.where(ftn_parents_all>=0)[0]

                if ftn_from not in graph_backward.keys():
                    graph_backward[ftn_from]=list(ftn_parents)
                else:
                    temp1=graph_backward.get(ftn_from)+list(ftn_parents)
                    graph_backward[ftn_from]=list(set(temp1))
        if new_fdg:
            self.graph_new=graph_backward
        else: self.graph_original=graph_backward

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
            The general method to get seqence from function_idx to leaf nodes in the graph
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
        seq_new = []
        path=[]
        visited = [False] * self.fdg.num_ftn
        _get_sequences(self.graph_original,ftn_idx,visited,path,seq_ori)

        # get sequences from new FDG (some functions may not appear in FDG)
        if ftn_idx in self.graph_new:
            path=[]
            visited = [False] * self.fdg.num_ftn
            _get_sequences(self.graph_new, ftn_idx, visited, path, seq_new)

        return {'original':seq_ori,'new':seq_new}

    def _get_shortest_sequences(self):
        ftn_shortest_seq_dict={}
        for ftn_idx in range(1,self.fdg.num_ftn):
            seqs={}
            sequences=self._generate_sequences(ftn_idx)
            # get the shortest sequence and length of shortest sequence
            for key,value in sequences.items():
                if len(value)>0:
                    ele_size=[len(item) for item in value]
                    min_size=min(ele_size)
                    shortest_seq=[item for item in value if len(item)==min_size]
                    seqs[key]={'sequences':shortest_seq, 'depth':min_size-1}

                else:
                    seqs[key]={'sequences':[], 'depth':0}

            ftn_shortest_seq_dict[ftn_idx]=seqs
        return ftn_shortest_seq_dict

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
                        seqs[key] = {'sequences': value, 'seq_size':len(value),'shortest':shortest_seq[0],'shortest_depth':min_size-1}

                    else:
                        seqs[key] = {'sequences': [], 'seq_size':0, 'shortest':[],'shortest_depth':0}

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

    def generate_target_sequences(self):
        seq_dict={}
        # get shortest sequences for all functions except constructor
        self._generate_backward_graph(True)  # build graph based on new fdg
        self._generate_backward_graph(False)  # build graph based on original fdg
        ftn_shortest_sequence_dict=self._get_shortest_sequences()
        print(f'shortest sequences:\n{ftn_shortest_sequence_dict}')
        for ftn_idx in self.ftn_idx_not_covered:
            l_p_dict=self._get_label_parent(ftn_idx)

            parenet_groups=[value for value in l_p_dict.values()]
            # case 1: only one parent or multiple parents connected by edges with same labels
            if len(parenet_groups)==1:
                for p_ftn in parenet_groups[0]:
                    p_sequences=[]
                    p_seq_dict=self._generate_sequences(p_ftn)
                    # only use sequences from new fdg, unreachable parent functions will never reach child functions
                    seqs=p_seq_dict['new']
                    seqs.sort(key=len,reverse=True) # sort based on size of sequences in descending
                    for seq in seqs:
                        p_sequences.append([(len(seq)-1,p_ftn),ftn_idx])

                    if ftn_idx not in seq_dict.keys():
                        seq_dict[ftn_idx]=p_sequences
                    else:
                        seq_dict[ftn_idx]+=p_sequences
            # case 2: multiple parents connected by edges with different labels
            else:
                # get sequences from single parent
                # consider the special case that target function can be reached by one parent directly
                parents=[]
                for values in l_p_dict.values():
                    parents+=values
                parents=list(set(parents))
                for p_idx in parents:
                    if ftn_shortest_sequence_dict[p_idx]['new']['depth']==self.fdg.depth_all_ftns_reached:

                        if ftn_idx not in seq_dict.keys():
                            seq_dict[ftn_idx]=[[(ftn_shortest_sequence_dict[p_idx]['new']['depth'],p_idx),ftn_idx]]
                        else:
                            seq_dict[ftn_idx]+=[[(ftn_shortest_sequence_dict[p_idx]['new']['depth'],p_idx),ftn_idx]]


                # get sequences through parent combination
                parent_combinations=self._get_combination(parenet_groups)
                collection_seq=[] # save the sequences generated

                for comb in parent_combinations:
                    # order elements in comb. sequences from new fdg with higher length has higher priority
                    comb_priority=[ftn_shortest_sequence_dict[ele]['new']['depth'] for ele in comb]

                    # if comb_priority==[0]*len(comb): continue
                    if 0 in comb_priority:continue  # if there is one sequence not reachable,ignore the combination
                    temp=np.array([comb,comb_priority]).T
                    temp=temp[temp[:,1].argsort()][::-1]
                    order_comb=temp[:,0]
                    order_comb_priority=temp[:,1]

                    # replace ele in order_comb with its shortest sequence
                    sequences=[]
                    for i in range(1,len(order_comb)):
                        if ftn_shortest_sequence_dict[order_comb[i]]['new']['sequences']:
                            p_seq=ftn_shortest_sequence_dict[order_comb[i]]['new']['sequences'][0]
                        else:
                            p_seq = ftn_shortest_sequence_dict[order_comb[i]]['original']['sequences'][0]
                        p_seq=p_seq[0:-1] # remove 0 from the sequence (the last element, denote constructor)
                        p_seq.reverse()
                        sequences.append(p_seq)

                    # merging sequence
                    if len(sequences)>=2:
                        temp=sequences[0]
                        for i in range(1,len(sequences)):
                            temp=merge_2_list(temp,sequences[i])
                    else:
                        temp=sequences[0]

                    # form final sequence
                    # the first two elements: the length of the first parent, parent index
                    final_seq=[(order_comb_priority[0],order_comb[0])]

                    # print(f'final_seq={final_seq}')
                    final_seq+=temp
                    final_seq.append(ftn_idx)

                    # print(f'final_seq={final_seq}')
                    # print(f'comb={comb}  size={len(comb)}')
                    collection_seq.append(final_seq)
                # limit the number of sequences
                if len(collection_seq)>self.num_seq_limit:
                    collection_seq.sort(key=len)
                    collection_seq=collection_seq[0:self.num_seq_limit]
                if ftn_idx in seq_dict.keys():
                    seq_dict[ftn_idx] += collection_seq
                else:
                    seq_dict[ftn_idx]=collection_seq
        print(f'sequences generated={seq_dict}')
        return seq_dict

    def generate_target_sequences_1(self):
        all_sequences={}
        # get shortest sequences for all functions except constructor
        self._generate_backward_graph(True)  # build graph based on new fdg
        self._generate_backward_graph(False)  # build graph based on origial fdg

        for ftn_idx in self.ftn_idx_not_covered:
            # get labels and parents
            l_p_dict = self._get_label_parent(ftn_idx)
            assert len(l_p_dict)>0

            #------------------------------------------------
            # case 1: consider each parent node individually
            # -----------------------------------------------
            parents = []
            for values in l_p_dict.values():
                parents += values
            parents = list(set(parents)) # remove repeated elements

            for p_idx in parents:
                p_sequences = []
                #get sequences for parent p_idx
                p_seq_dict = self._generate_sequences(p_idx)

                # only use sequences from new fdg,unreachable parent functions will never reach child functions
                seqs = p_seq_dict['new']
                # seqs.sort(key=len, reverse=True)  # sort based on sequence length in descending
                for seq in seqs:
                    # only consider sequences with depth=self.fdg.depth_all_ftns_reached
                    # len of a sequence is 1 more than depth of the sequence
                    # due to constructor included in sequence
                    if len(seq) - 1 == self.fdg.depth_all_ftns_reached:
                        p_sequences.append([([self.fdg.depth_all_ftns_reached], p_idx), ftn_idx])

                if ftn_idx not in all_sequences.keys():
                    all_sequences[ftn_idx] = p_sequences
                else:
                    all_sequences[ftn_idx] += p_sequences

            # ------------------------------------------------
            # case 2: consider combined parent nodes
            # ------------------------------------------------
            parenet_groups = [value for value in l_p_dict.values()]
            if len(parenet_groups)>1:
                # do combinations

                # key-value format: target function: [[target function, other function, constructor],...]
                ftn_seq_and_shortest_seq_dict = self._get_sequences_and_shortest_sequences()

                # get sequences through parent combination
                parent_combinations=self._get_combination(parenet_groups)

                collection_seq=[] # save the sequences generated
                for comb in parent_combinations:
                    # order elements in comb.
                    # sequences from new fdg with higher length has higher priority

                    comb_priority_shortest_seq_len=[ftn_seq_and_shortest_seq_dict[ele]['new']['shortest_depth'] for ele in comb]
                    comb_priority_seq_size=[ftn_seq_and_shortest_seq_dict[ele]['new']['seq_size'] for ele in comb]

                    # if comb_priority==[0]*len(comb): continue
                    # in this version, unreachable sequences are not considered
                    if 0 in comb_priority_seq_size:continue  # if there is one sequence not reachable,ignore the combination
                    temp=np.array([comb,comb_priority_seq_size,comb_priority_shortest_seq_len]).T
                    # temp=temp[temp[:,1].argsort()][::-1]
                    temp = temp[np.lexsort((temp[:, 1],temp[:,2]))][::-1]
                    order_comb=temp[:,0]



                    # replace ele in order_comb with its shortest sequence
                    # the first ele is not replaced as it is used to find base states and already executed.
                    sequences=[]
                    for i in range(1,len(order_comb)):
                        if ftn_seq_and_shortest_seq_dict[order_comb[i]]['new']['shortest']:
                            p_seq=ftn_seq_and_shortest_seq_dict[order_comb[i]]['new']['shortest']

                            # remove 0 from the sequence (the last element denotes constructor)
                            p_seq=p_seq[0:-1]
                            p_seq.reverse()
                            sequences.append(p_seq)

                    if len(sequences)==0:
                        continue
                    # merging sequence
                    if len(sequences)>=2:
                        temp=sequences[0]
                        for i in range(1,len(sequences)):
                            temp=merge_2_list(temp,sequences[i])
                    else:
                        temp=sequences[0]


                    # form final sequence: first tuple+ rest sequence+ target function
                    # handle the first ele in comb
                    seq_1st_ele=ftn_seq_and_shortest_seq_dict[order_comb[0]]['new']['sequences']

                    len_ary=[len(seq)-1 for seq in seq_1st_ele]
                    len_ary=list(set(len_ary))
                    len_ary.sort(reverse=True) #order in descending order
                    len_ary=[ele for ele in len_ary if ele<=self.fdg.depth_all_ftns_reached]
                    if len(len_ary)>=2:
                        # format: ([depths that parent ftn can be reached],parent_ftn_idx)
                        final_seq = [(len_ary[0:2], order_comb[0])]
                    else:
                        final_seq = [([len_ary[0]], order_comb[0])]


                    # print(f'final_seq={final_seq}')
                    final_seq+=temp
                    final_seq.append(ftn_idx)

                    collection_seq.append(final_seq)

                # after going through all combinations
                # limit the number of sequences
                if len(collection_seq)>self.num_seq_limit:
                    collection_seq.sort(key=len) # in ascending order
                    collection_seq=collection_seq[0:self.num_seq_limit]
                if ftn_idx in all_sequences.keys():
                    all_sequences[ftn_idx] += collection_seq
                else:
                    all_sequences[ftn_idx]=collection_seq
        print(f'sequences generated 1={all_sequences}')
        return all_sequences



    def info_1st_and_following_sequence_execution(self,sequnces_dict:dict):
        info_1st_step={}
        for key,sequences in sequnces_dict.items():
            if not sequences:continue
            for seq in sequences:
                ele_1st=seq[0] # format: ([depth])
                ele_2nd=seq[1] # the first ftn to be executed
                for depth in ele_1st[0]:
                    if depth not in info_1st_step.keys():
                        info_1st_step[depth]={ele_1st[1]:[ele_2nd]}
                    else:
                        if ele_1st[1] not in info_1st_step[depth].keys():
                            info_1st_step[depth][ele_1st[1]]=[ele_2nd]
                        else:
                            info_1st_step[depth][ele_1st[1]] += [ele_2nd]

        info_following_step = {}
        for key,sequences in sequnces_dict.items():
            if not sequences: continue
            for seq in sequences:
                if len(seq)>2:
                    for i in range(2,len(seq)):
                        pre_ftn=seq[i-1]
                        ftn=seq[i]
                        if i not in info_following_step.keys():
                            info_following_step[i]={pre_ftn:[ftn]}
                        else:
                            if pre_ftn not in info_following_step[i].keys():
                                info_following_step[i][pre_ftn]=[ftn]
                            else:info_following_step[i][pre_ftn]+=[ftn]
        return info_1st_step,info_following_step


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

    def get_sequences_by_level(self,level:int):
        """

        :param level: number of parent groups considered
        :return:
        """
        assert self.ftn_seq_and_shortest_seq_dict
        all_sequences={}
        if level<1: return  all_sequences
        if level==1:
            # ------------------------------------------------
            # levl 1: consider only one parent node
            # -----------------------------------------------
            for ftn_idx in self.ftn_idx_not_covered:
                # get labels and parents
                l_p_dict = self._get_label_parent(ftn_idx)
                if len(l_p_dict) == 0: continue

                parents = []
                for values in l_p_dict.values():
                    parents += values
                parents = list(set(parents))  # remove repeated elements

                for p_idx in parents:
                    p_sequences = []
                    # get sequences for parent p_idx
                    p_seq_dict = self._generate_sequences(p_idx)

                    # only use sequences from new fdg,unreachable parent functions will never reach child functions
                    seqs = p_seq_dict['new']
                    # seqs.sort(key=len, reverse=True)  # sort based on sequence length in descending
                    for seq in seqs:
                        # only consider sequences with depth=self.fdg.depth_all_ftns_reached
                        # len of a sequence is 1 more than depth of the sequence
                        # due to constructor included in sequence
                        if len(seq) - 1 == self.fdg.depth_all_ftns_reached:
                            p_sequences.append([([self.fdg.depth_all_ftns_reached], p_idx), ftn_idx])

                    if ftn_idx not in all_sequences.keys():
                        all_sequences[ftn_idx] = p_sequences
                    else:
                        all_sequences[ftn_idx] += p_sequences

        else:
            # ------------------------------------------------
            # level 2 and higher : consider combined parent nodes
            # ------------------------------------------------
            for ftn_idx in self.ftn_idx_not_covered:
                # get labels and parents
                l_p_dict = self._get_label_parent(ftn_idx)
                if len(l_p_dict) == 0: continue

                parent_groups = [values for values in l_p_dict.values()]
                if len(parent_groups)>=level:
                    # get sequences through parent combination
                    parent_combinations = self._get_combination(parent_groups, level)

                    collection_seq = []  # save the sequences generated
                    for comb in parent_combinations:
                        # order elements in comb.
                        # sequences from new fdg with higher length has higher priority

                        comb_priority_shortest_seq_len = [
                            self.ftn_seq_and_shortest_seq_dict[ele]['new']['shortest_depth'] for ele in comb]

                        # in this version, unreachable sequences are not considered
                        if 0 in comb_priority_shortest_seq_len: continue  # if there is one sequence not reachable,ignore the combination
                        temp = np.array([comb, comb_priority_shortest_seq_len]).T
                        # temp=temp[temp[:,1].argsort()][::-1]
                        temp = temp[temp[:, 1].argsort()][::-1]
                        # temp = temp[np.lexsort((temp[:, 1]))][::-1]
                        order_comb = temp[:, 0]

                        # replace elements in order_comb with its shortest sequence
                        sequence=[]
                        for ele in order_comb:
                            seq=self.ftn_seq_and_shortest_seq_dict[ele]['new']['shortest'][0:-1]
                            seq.reverse()
                            sequence.append(seq)

                        p_first = sequence[0]
                        p_first_length = len(sequence[0])
                        merge_seq = []
                        if len(sequence) == 2:
                            merge_seq = merge_fix_list_1(sequence[0], sequence[1])
                        elif len(sequence) > 2:
                            if p_first_length == 1:
                                for i in range(2, len(sequence)):
                                    merge_seq = merge_2_list(merge_seq, sequence[i])
                            else:
                                for i in range(2, len(sequence)):
                                    merge_seq = merge_fix_list_1_specified_lenth(merge_seq, sequence[i],
                                                                                 p_first_length)

                        # form final sequence
                        final_seq_2nd_part = merge_seq[p_first_length:]  # remove the elements from the first parent
                        final_seq_2nd_part.append(ftn_idx)
                        final_seq = [([p_first_length], order_comb[0])]+ final_seq_2nd_part

                        collection_seq.append(final_seq)

                    # after going through all combinations
                    # limit the number of sequences
                    if len(collection_seq) > self.num_seq_limit:
                        collection_seq.sort(key=len)  # in ascending order
                        collection_seq = collection_seq[0:self.num_seq_limit]

                    if ftn_idx in all_sequences.keys():
                        all_sequences[ftn_idx] += collection_seq
                    else:
                        all_sequences[ftn_idx] = collection_seq

        return all_sequences

    def organize_sequences_for_execution_by_level(self,sequences_dict,level:int):

        info_1st_step = {}
        info_following_step = {}
        if level<1: return info_1st_step,info_following_step

        for key, sequences in sequences_dict.items():
            if not sequences: continue
            for seq in sequences: # [([depths],p_ftn),ftn1,ftn2,,,]
                ele_1st = seq[0]  # format: ([depth])
                ele_2nd = seq[1]  # the first ftn to be executed
                for depth in ele_1st[0]:
                    if depth not in info_1st_step.keys():
                        info_1st_step[depth] = {ele_1st[1]: [ele_2nd]}
                    else:
                        if ele_1st[1] not in info_1st_step[depth].keys():
                            info_1st_step[depth][ele_1st[1]] = [ele_2nd]
                        else:
                            if ele_2nd not in info_1st_step[depth][ele_1st[1]]:
                                info_1st_step[depth][ele_1st[1]] += [ele_2nd]

        if level>1:
            for key, sequences in sequences_dict.items():
                if not sequences: continue
                for seq in sequences:
                    if len(seq) > 2:
                        for i in range(2, len(seq)):
                            pre_ftn = seq[i - 1]
                            ftn = seq[i]
                            if i not in info_following_step.keys():
                                info_following_step[i] = {pre_ftn: [ftn]}
                            else:
                                if pre_ftn not in info_following_step[i].keys():
                                    info_following_step[i][pre_ftn] = [ftn]
                                else:
                                    if ftn not in info_following_step[i][pre_ftn]:
                                        info_following_step[i][pre_ftn] += [ftn]
        return info_1st_step, info_following_step


# merge two lists, common elements are considered once
# the order in each list is maintained in combined list
# the first element in the combined result should be the first element in a
def merge_2_list(a,b)->list:
    a1 = [] # temporarily save elements in a that are not matched with element in b
    b1 = [] # temporarily save elements in b that are not matched with element in a
    c = [] # save the combined result

    if not b:
        c=a
        return c
    if not a:
        c=b
        return c
    j = 0
    while j < len(b):
        flag = True
        for i in range(len(a)):
            if b[j] == a[i]:
                # add elements to c
                if len(a1) > 0:
                    c += a1
                    a1 = []
                if len(b1) > 0:
                    c += b1
                    b1 = []
                c.append(a[i])
                # update a
                a = a[i + 1:]  # remove elements that already moved to c
                # check if either a or b reaches end
                if len(a) == 0:
                    if j < len(b) - 1:
                        c += b[j + 1:]
                    return c
                else:
                    if j == len(b) - 1:
                        c += a
                        return c
                    else: # both a and b still do not reach end, go to the next iteration
                        j += 1
                flag = False
                break
            else: # do not match, save a[i] to a1
                a1.append(a[i])

        if flag:
            if j < len(b) - 1:
                a1 = []
                b1.append(b[j])
                j += 1
            else:
                c+=a
                c.append(b[j])
                return c

def merge_fix_list_1(list_1,list_2):
    c=list_1
    j=0
    while j < len(list_2):
        for i in range(len(list_1)):
            if j==len(list_2):break
            if list_1[i]==list_2[j]:
                j+=1

        if j<len(list_2):
            c+=list_2[j:]
        break
    return c


def merge_fix_list_1_specified_lenth(list_1,list_2,specified_length):
    c=list_1[0:specified_length]
    j=0
    while j < len(list_2):
        for i in range(specified_length):
            # if i==len(list_2)-1
            if list_1[i]==list_2[j]:
                j+=1
        break
    a=list_1[specified_length-1:]
    b=list_2[j:]
    re=merge_2_list(a,b)
    if re:
        c+=re[1:]
    return c








if __name__=='__main__':
    # # ftn_info=Function_info('/home/wei/PycharmProjects/Contracts/_wei/Crowdsale.sol', 'Crowdsale')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    #
    # function_dict=ftn_info.functions_dict_slither()
    #
    # fdg = FDG(function_dict)
    # function_list = [1, 5, 7, 6, 10]
    # fdg.build_fdg_3d_array(function_list)





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

    seq= Sequence()
    a=[1,2,5,7]
    b=[4,5]
    print(merge_2_list(a,b))

    print(merge_fix_list_1([1,2,5,7],[4,5]))
    print(merge_fix_list_1([1, 2, 5, 7], [2,8]))
    print(merge_fix_list_1([1, 2, 5, 7], [2,7]))
    print(merge_fix_list_1([1, 2, 5, 7], [2,1]))
    print(merge_fix_list_1([1, 2, 5, 7], [1,3,7]))
    print(merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1, 5],3))










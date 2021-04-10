import numpy as np
import itertools as it

from fdg.FDG import FDG
s=[]

class Sequence():
    def __init__(self, fdg=None,target_ftn=None):
        self.fdg=fdg
        self.ftn_name=target_ftn # target function
        self.depths_have_multi_in_edges=[] # save the depth values for self.ftn_name at which there are multiple incoming edges
        self.parents_at_depth_has_multi_in_edges=[] # the nodes(functions) at each depth with multiple incoming edges

    def get_combination(self, list_for_comb):
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

    def get_parents_combination(self)->dict:
        """
        based on original fdg, which is not updated
        :return:
        """
        combination_dict={}
        for depth in range(1,self.fdg.depth_all_ftns_reached):
            ftn_index=self.fdg.ftn_to_index[self.ftn_name]
            ftn_parents=self.fdg.original_matrix_fdg[depth,ftn_index,:]
            unique_labels=np.unique(ftn_parents[ftn_parents>=0])

            # do combination when have two or more edges with different labels
            if len(unique_labels)>=2:
                # group based on the label values
                groups=[list(np.where(ftn_parents== lbl)[0]) for lbl in unique_labels]
                combination_dict[depth]=self.get_combination(groups)
                self.depths_have_multi_in_edges.append(depth)
                self.parents_at_depth_has_multi_in_edges.append(list(np.where(ftn_parents >= 0)[0]))

        return combination_dict


    def build_graph_backward(self, depth, function_idx):
        """
        -find a graph dict for function_idx to depth 1
        -breadth first search
        -There is a case, that the target function can appear at lower depth. So in this case, just ignore it
        :param depth: the starting depth
        :param function_idx:  the starting function's index
        :return:
        """
        graph = {}
        ftns_list = [function_idx]
        for d in range(depth, 0, -1):
            new_ftns = set()
            for ftn_idx in ftns_list:
                ftns_parents = list(np.where(self.fdg.original_matrix_fdg[d, ftn_idx, :] >= 0)[0])
                # remove the target function if it containes
                if function_idx in ftns_parents:
                    ftns_parents.remove(function_idx)
                graph[ftn_idx] = ftns_parents
                new_ftns = new_ftns.union(set(ftns_parents))
            ftns_list = list(new_ftns)  # prepare for the next loop
        return graph

    def get_sequences_general(self, graph, function_idx, visited, path, se):
        """
        The general method to get seqence from function_idx to leaf nodes in the graph
        :param graph: a dict
        :param function_idx: the source
        :param visited:
        :param path:
        :param se: save the resulted sequences
        :return:
        """
        # checking all the visited nodes
        visited[function_idx] = True
        path.append(function_idx)

        if function_idx not in graph:  # leaf nodes
            # se.append(list(np.copy(path))[1::][::-1])
            se.append(list(np.copy(path)))
        else:
            for i in graph[function_idx]:
                if visited[i] == False:
                    self.get_sequences_general(graph, i, visited, path, se)
        path.pop()
        visited[function_idx] = False

    def dict_ftns_short_sequences(self):
        function_idx = self.fdg.ftn_to_index[self.ftn_name]

        # build the graph for target function
        dict_shortest = {}
        # for i in range(1, self.fdg.matrix_fdg.shape[0]):
        for i in range(1, self.fdg.original_matrix_fdg.shape[0]):
            graph = self.build_graph_backward(i, function_idx)
            if len(graph) > 0:
                # generte sequences
                visited = [False] * self.fdg.num_ftn
                path = []
                se = []
                seq = self.get_sequences_general(graph, function_idx, visited, path, se)
                if len(se) > 0:
                    dict_i = {}
                    for s in se:
                        s = s[1::][::-1]
                        if s[-1] in dict_i.keys():
                            dict_i[s[-1]] += [s]
                        else:
                            dict_i[s[-1]] = [s]
                    dict_shortest[i] = dict_i

        return dict_shortest




    def flat_2d_list_with_2_elements(self,list_2d_2element):
        """
               :param list_2d: [[2,6],[2,7,8]]
               :return: [2,6,7,8]
               """
        re_list = []
        c_len = [len(c) for c in list_2d_2element]
        # convert the 2d list to 2d matrix with the second dimension length to min length of its elements
        np_a = np.array([list_2d_2element[0][0:min(c_len)]])
        np_a = np.concatenate((np_a, np.array([list_2d_2element[1][0:min(c_len)]])), axis=0)

        # compare by element-wise substract
        b = np.subtract(np_a, np_a[0, :])
        for j in range(b.shape[1]):
            flag = True
            if len(np.where(b[:, j] > 0)[0]) > 0:
                flag = False
                re_list = list_2d_2element[0]
                if len(list_2d_2element[1])>j+1:
                    re_list += list_2d_2element[1][j + 1:len(list_2d_2element[1])]
                break

        if flag:  # the case all elements share the prefix with length min(c_len)
            re_list = list_2d_2element[0]
            if len(list_2d_2element[1]) > min(c_len) + 1:
                re_list += list_2d_2element[1][j + 1:len(list_2d_2element[1])]

        return re_list

    def flat_2d_list_with_more_elements(self, list_2d_more_elements):
        result=[]
        len_list=len(list_2d_more_elements)
        for i in range(len_list-1):
            list_2d_2e=[list_2d_more_elements[i],list_2d_more_elements[i+1]]
            re=self.flat_2d_list_with_2_elements(list_2d_2e)
            result+=re
        return result

    def generate_depth_seqeunces(self):
        """
        :param combs: {depth: [combinations] }
        :param ftns_seq_short:
        :return:
        """
        seq_dict = {}
        combs = self.get_parents_combination()
        ftns_seq_short = self.dict_ftns_short_sequences()

        for depth in self.depths_have_multi_in_edges:
            parent_combs = combs[depth]
            s_collect=[]
            for p_comb in parent_combs:
                starting_ftn = p_comb[0]
                # find the starting function in the combination
                for p_ftn in p_comb:
                    # if p_ftn only appeara at the current depth, then it can be the starting function
                    flag_start = True
                    for d in range(1, depth):
                        if p_ftn in ftns_seq_short[d].keys():
                            flag_start = False
                            break
                    if flag_start:
                        starting_ftn = p_ftn
                        break

                other_ftns = []
                if len(p_comb) >= 3:  # need to consider the prefix among functions except the first one
                    # replace p_ftn with its shortest sequence to the target functio
                    for p_ftn in p_comb:
                        if p_ftn != starting_ftn:
                            for j in range(1, depth + 1):
                                if p_ftn in ftns_seq_short[j].keys():
                                    other_ftns.append(ftns_seq_short[j][p_ftn])
                                    break  # only take the shortest one (that means once found from lower depth, then stop)
                    # handle prefix among other functions execpt for the starting function
                    # example of  other_ftns:[[7,9],[7,10]]
                    other_ftns=self.flat_2d_list_with_more_elements(other_ftns)
                else:
                    other_ftns = p_comb[1::]

                # ----------------------------------------
                # assume other_ftns is one-dimension list
                # remove prefix between the sequence related to starting_ftn and other_ftns
                seq_of_starting_ftn = ftns_seq_short[depth][starting_ftn][0]  # To-Do
                min_len = min(len(other_ftns), len(seq_of_starting_ftn))
                for i in range(min_len):
                    if other_ftns[i] != seq_of_starting_ftn[i]:
                        other_ftns = other_ftns[i::]
                        break

                # form sequences
                seq = [starting_ftn]
                seq += other_ftns
                seq += [self.fdg.ftn_to_index[self.ftn_name]]
                s_collect.append(seq)
            seq_dict[depth]=s_collect

        return seq_dict

    def final_depth_sequences(self):
        combs = self.get_parents_combination()

        ftns_seq_short = self.dict_ftns_short_sequences()
        # print(f'ftns and shortest sequences:{ftns_seq_short}')
        result = self.generate_depth_seqeunces(combs, ftns_seq_short)
        return result

if __name__=='__main__':
    functions_dict={'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'], 'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'], 'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'], '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'], 'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'], 'f6': ['decreaseApproval', [], ['allowed'], '66188463'], 'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'], 'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'], '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'], 'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'], 'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'], 'f0': ['constructor', [], ['owner'], '8afc3605']}

    fdg=FDG(functions_dict)
    seq=Sequence(fdg,'burn')
    # re=seq.generate_depth_seqeunces()
    # re=seq.dict_ftns_short_sequences()
    re=seq.get_parents_combination()
    # dict_shortest=seq.dict_ftns_short_sequences()
    # print(dict_shortest)
    # seq.generate_depth_seqeunces(combs,dict_shortest)
    # # re= seq.get_sequences_at_depth(2, "burn")
    print(re)





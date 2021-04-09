import numpy as np
import itertools as it

from fdg.FDG import FDG
s=[]

class Sequence():
    def __init__(self, fdg=None,target_ftn=None):
        self.fdg=fdg
        self.ftn_name=target_ftn # target function
        self.depths_multi_in_edges=[] # save the depth values for self.ftn_name at which there are multiple incoming edges
        self.ftns_at_multi_in_edges=[] # the nodes(functions) at each depth with multiple incoming edges


    def get_combination_sequences_for_ftn(self)->dict:
        """
        based on original fdg, which is not updated
        :return:
        """
        combination_dict={}
        for depth in range(1,self.fdg.depth_all_ftns_reached+1):
            ftn_index=self.fdg.ftn_to_index[self.ftn_name]
            ftn_dependency_status=self.fdg.original_matrix_fdg[depth,ftn_index,:]
            unique_labels=np.unique(ftn_dependency_status[ftn_dependency_status>=0])

            if len(unique_labels)>=2: # do combination
                # group based on the label values
                groups=[list(np.where(ftn_dependency_status== lbl)[0]) for lbl in unique_labels]
                combination_dict[depth]=self.get_combination_sequences(groups)
                self.depths_multi_in_edges.append(depth)
                self.ftns_at_multi_in_edges.append(list(np.where(ftn_dependency_status>=0)[0]))

        return combination_dict

    def get_combination_sequences(self,list_for_comb):
        """
        group indices with the same value [[1,4], [2,6],[5]]
        from each group, choose 1 index
        generate combinations with length from 2 to the number of groups
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



    def get_shortest_sequences(self)->dict:
        """
        based on updated fdg
        generate shortest sequences for some functions to reach the target function
        :return:
        """
        dict_re = {}
        if len(self.depths_multi_in_edges)>0 and len(self.ftns_at_multi_in_edges)>0:
            # get all possible sequences for eqch function at possible depths
            for i in range(len(self.depths_multi_in_edges)):
                sequences=self.get_sequences_for_ftn(self.depths_multi_in_edges[i],self.ftn_name)
                dict_re[self.depths_multi_in_edges[i]]=sequences

            # heuristically shorten sequences at depth 1
            ftn_content_at_depth_1 = self.fdg.original_matrix_fdg[1, self.fdg.ftn_to_index[self.ftn_name], :]
            #all functions that reaches target function at the depth 1
            ftns=np.where(ftn_content_at_depth_1>=0)[0]
            for key,value in dict_re.items():
                if isinstance(value,dict):
                    for k,v in value.items():
                        if k in ftns and len(v)>0:
                            dict_re[key][k]=[k]
        return dict_re

    def get_sequences_for_ftn(self,depth,ftn_name)->dict:
        ftn_index=self.fdg.ftn_to_index[ftn_name]
        visited=[False]*self.fdg.num_ftn
        path=[]
        se=[]
        # build the graph
        graph=self.build_graph(depth,ftn_index)

        # generate sequence
        self.get_sequences(graph,ftn_index,visited,path,se)

        # convert sequence to dict
        re_dict={}
        for s in se:
            if s[-1] in re_dict:
                re_dict[s[-1]]+=[s]
            else:re_dict[s[-1]]=[s]
        return re_dict

    def get_sequences(self, graph,ftn_idx,visited, path,se):
        # checking all the visited nodes
        visited[ftn_idx] = True
        path.append(ftn_idx)

        if ftn_idx not in graph:
            #print(path)
            se.append(list(np.copy(path))[1::][::-1])

        else:
            for i in graph[ftn_idx]:
                if visited[i] == False:
                   self.get_sequences(graph, i, visited, path,se)
        path.pop()
        visited[ftn_idx] = False

    def build_graph(self,depth,ftn_idx):
        graph={}
        ftns_list=[ftn_idx]
        for d in range(depth,0,-1):
            new_ftns=set()
            for idx in ftns_list:
                ftns_depedent=list(np.where(self.fdg.matrix_fdg[d,idx,:]>=0)[0])
                graph[idx]=ftns_depedent
                new_ftns=new_ftns.union(set(ftns_depedent))
            ftns_list=list(new_ftns) # prepare for the next loop
        return graph



    def generate_depth_seqeunces_for_ftn(self,combs,ftns_seq_short):
        """
        :param combs: {depth: [combinations] }
        :param ftns_seq_short:
        :return:
        """
        seq_dict={}
        for depth in self.depths_multi_in_edges:
            s_collect=[]
            comb_=combs[depth]
            ftns_=ftns_seq_short[depth]

            for item in comb_:
                # ========================================
                # replace each element in the item with its shortest sequence to the target function
                s = []
                item_len = []
                for ftn_i in item:
                    if ftn_i not in ftns_:
                        item_len.append(0)
                    else:
                        item_len.append(len(ftns_[ftn_i]))
                if 0 in item_len:  # means that there is one function not reachable
                    continue

                s=[[item[0]]]
                for i in range(1, len(item)):
                    s.append(ftns_[item[i]])
                comb_re = [com for com in it.product(*s)]
                # s=[]
                # item_len=[]
                # for ftn_i in item:
                #     if ftn_i not in ftns_:
                #         item_len.append(0)
                #     else:item_len.append(len(ftns_[ftn_i]))
                #
                # if 0 in item_len: # means that there is one function not reachable
                #     continue
                #
                # if item_len[0]==max(item_len):
                #     s=[[item[0]]] # starting ftn
                #     # replace the left functions with their shortest sequence to the target function
                #     for i in range(1,len(item)):
                #         s.append(ftns_[item[i]])
                # else:
                #     max_idx=item_len.index(max(item_len))
                #     s=[[item[max_idx]]] # starting ftn
                #     # replace the left functions with their shortest sequence to the target function
                #     for i in range(len(item)):
                #         if i!=max_idx:
                #             s.append(ftns_[item[i]])
                # # s.append([self.fdg.ftn_to_index[self.ftn_name]])  # add target function
                # comb_re=[com for com in it.product(*s)]

                #========================================
                # combined result handle
                # consider the shared prefix among multiple lists (the first several elements of these lists have the same values
                comb_re=list(comb_re[0])
                comb_t = [comb_re[0]]  # related the first element in the item(combination)
                sub_s = [] # related to other elements except the first one in the item
                if len(comb_re)>=3:
                    c_len = [len(c) for c in comb_re[1::]]

                    np_a = np.array([comb_re[1][0:min(c_len)]])
                    for i in range(2,len(comb_re)):
                        np_a=np.concatenate((np_a,np.array([comb_re[i][0:min(c_len)]])),axis=0)
                    b=np.subtract(np_a,np_a[0,:])

                    for j in range(b.shape[1]):
                        flag=True
                        if len(np.where(b[:,j]>0)[0])>0:
                            flag=False
                            sub_s=comb_re[1]
                            for m in range(2,len(comb_re)):
                                sub_s+=comb_re[m][j+1:len(comb_re[m])]
                            break

                    if flag:
                        sub_s = comb_re[1]
                        for m in range(2, len(comb_re)):
                            sub_s += comb_re[m][min(c_len):len(comb_re[m])]
                else:
                    if isinstance(comb_re[1],list):
                        sub_s +=comb_re[1]
                    else:sub_s +=[comb_re[1]]

                # handle the shared ftn elements related to the first element in the item
                # with the  ftn elements of the combined elements of the item
                min_len=min(len(sub_s),len(ftns_[comb_t[0]][0]))
                for i in range(min_len):
                    if sub_s[i]!=ftns_[comb_t[0]][0][i]:
                        sub_s=sub_s[i::]
                        break

                # form sequence
                comb_t+=sub_s
                comb_t+=[self.fdg.ftn_to_index[self.ftn_name]]  # add the target function to form the complete sequence
                s_collect.append(comb_t)

            seq_dict[depth]=s_collect
        return seq_dict



    def final_depth_sequences(self):
        combs = self.get_combination_sequences_for_ftn()
        print(f'combinations:{combs}')
        ftns_seq_short = self.get_shortest_sequences()
        print(f'ftns and shortest sequences:{ftns_seq_short}')
        result = self.generate_depth_seqeunces_for_ftn(combs, ftns_seq_short)
        return result

if __name__=='__main__':
    functions_dict={'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'], 'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'], 'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'], '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'], 'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'], 'f6': ['decreaseApproval', [], ['allowed'], '66188463'], 'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'], 'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'], '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'], 'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'], 'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'], 'f0': ['constructor', [], ['owner'], '8afc3605']}

    fdg=FDG(functions_dict)
    # seq=Sequence(fdg,'burn')
    # combs = seq.get_combination_sequences_for_ftn()
    # # do some updation
    # seq.fdg.matrix_fdg[1, 3, 5] = -3
    # seq.fdg.matrix_fdg[1, 3, 6] = -3
    # seq.fdg.matrix_fdg[1, 11, 10] = -3
    # ftns_seq_short = seq.get_shortest_sequences()
    # print(seq.generate_depth_seqeunces_for_ftn(combs,ftns_seq_short))

    seq = Sequence(fdg, 'burn')
    combs = seq.get_combination_sequences_for_ftn()
    # do some updation
    seq.fdg.matrix_fdg[1, 3, 5] = -3
    seq.fdg.matrix_fdg[1, 3, 6] = -3
    seq.fdg.matrix_fdg[1, 11, 10] = -3
    seq.fdg.matrix_fdg[2, 2, 11] = -3
    seq.fdg.matrix_fdg[2, 3, 11] = -3
    seq.fdg.matrix_fdg[2, 8, 11] = -3
    seq.fdg.matrix_fdg[2, 2, 3] = -3
    seq.fdg.matrix_fdg[2, 8, 3] = -3
    seq.fdg.matrix_fdg[2, 11, 3] = -3
    seq.fdg.matrix_fdg[2, 11, 10] = -3
    seq.fdg.matrix_fdg[2, 11, 8] = -3
    seq.fdg.matrix_fdg[2, 3, 8] = -3
    seq.fdg.matrix_fdg[2, 3, 9] = -3

    # re=seq.get_sequences_for_ftn(3,'transferFrom')
    # ftns_seq_short = seq.get_shortest_sequences()
    # print(seq.generate_depth_seqeunces_for_ftn(combs, ftns_seq_short))


    # # re=seq.get_shortest_sequences()
    re=seq.final_depth_sequences()
    print(re)








from copy import copy
import numpy as np
import csv
import sys
import ast

from fdg import FDG_global
from fdg.funtion_info import Function_info
import json

class FDG():

    def __init__(self,functions_dict):
        """
        :param functions_dict:
        key: [function_name,read_sv_list,write_sv_list,identifier]
         {'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b',1], ...}
        """
        self.limit_maximum_depth=5
        self.ftn_to_index={}
        self.index_to_ftn={}
        self.ftn_0_to_index={} # not function full name ( i.e., without parentheses and argument types)
        self.index_to_ftn_0={}
        self.ftn_to_selector={}
        self.index_to_selector={}
        self.labels=set()
        self.label_to_index={}
        self.index_to_label={}
        self.num_ftn=0
        self.num_label=0

        for key, value in functions_dict.items():
            #value: full_name,read_set, write_set,hash_value, index
            self.ftn_to_index[value[0]]=value[4]
            self.ftn_0_to_index[str(value[0]).split('(')[0]]=value[4]
            self.index_to_ftn[value[4]]=value[0]
            self.index_to_ftn_0[value[4]] = str(value[0]).split('(')[0]
            self.index_to_selector[value[4]] = value[3]
            self.ftn_to_selector[value[0]] = value[3]
            self.labels=self.labels.union(value[1]+value[2])

        self.labels = sorted(self.labels)
        for idx,item in enumerate(self.labels):
            self.index_to_label[idx]=item
            self.label_to_index[item]=idx

        self.num_ftn=len(list(functions_dict.keys()))
        self.num_label=len(self.labels)


        # create read, write matrices for all functions
        self.sv_read =-np.ones((self.num_ftn,self.num_label),dtype=int)
        self.sv_write =-np.ones((self.num_ftn,self.num_label),dtype=int)
        for key,value in functions_dict.items():
            read=-np.ones(self.num_label,dtype=int)
            write = -np.ones(self.num_label,dtype=int)

            if len(value[1])>0:
                for label in value[1]:
                    lbl_index=self.label_to_index[label]
                    read[lbl_index] = int(lbl_index)

            if len(value[2]):
                for label in value[2]:
                    lbl_index = self.label_to_index[label]
                    write[lbl_index] = int(lbl_index)

            if value[4]==0: #0: index of constructor()
                self.sv_read[0,:]=-np.ones(self.num_label,dtype=int).T
                self.sv_write[0,:] = -np.ones(self.num_label,dtype=int).T
            else:
                self.sv_read[value[4],:]=read.T
                self.sv_write[value[4],:] = write.T

        # get all edges
        self.nodes = list(range(self.num_ftn))
        self.edges={}
        self.nodes_wo_DD_edges=[]
        self.nodes_w_DD_edges=[]
        self.start_nodes=[]

        self.graph = {}
        self.graph[0]=list(range(1,self.num_ftn))
        # save edges
        for i in range(1,self.num_ftn):
            self.edges['f'+str(0)+",f"+str(i)]='none' # no edge label

        # ignore constructor
        for ftn_from in range(1, self.num_ftn):
            sv_write = self.sv_write[ftn_from, :]
            sv_write_idx = sv_write[sv_write >= 0]  # indices of state variables written
            if len(sv_write_idx)==0:
                continue
            for sv_w_idx in sv_write_idx:
                ftn_to = self.sv_read[:, sv_w_idx]
                ftn_to_idx = np.where(ftn_to >= 0)[0]  # indices of functions reading sv_w_idx
                if len(ftn_to_idx)==0:
                    continue
                for ftn_to in ftn_to_idx:
                    if ftn_from != ftn_to:
                        if ftn_from in self.graph.keys():
                            self.graph[ftn_from] += [ftn_to]
                        else:
                            self.graph[ftn_from] = [ftn_to]

                        if ftn_to not in self.nodes_w_DD_edges:
                            self.nodes_w_DD_edges.append(ftn_to)

                        # save edges
                        if 'f'+str(ftn_from)+",f"+str(ftn_to) not in self.edges.keys():
                            self.edges['f'+str(ftn_from)+",f"+str(ftn_to)]=self.index_to_label[sv_w_idx]


                # save nodes that have edges
                if ftn_from not in self.nodes_w_DD_edges:
                    self.nodes_w_DD_edges.append(ftn_from)

                # remove repeated element
                if ftn_from in self.graph.keys():
                    self.graph[ftn_from]=list(set(self.graph[ftn_from]))

        self.nodes_wo_DD_edges = list(set(self.nodes).difference(set(self.nodes_w_DD_edges)))
        self.nodes_wo_DD_edges=[item for item in self.nodes_wo_DD_edges if item not in [0,1]]








    def build_fdg_3d_array(self,function_list:list):
        """
        calculate the depth all functions are reached at least once.
        stop build FDG when all edges are included
        :param function_list: [10,3,6]
        :return:
        """
        # assert (len(function_list) > 0)
        edge_not_reached = copy(self.edge_list)
        ftn_mark = [False] * self.num_ftn
        ftn_mark[0] = True  # for constructor,set it true
        self.nodes_w_edges.append(0)
        self.graph[0] = function_list

        for ftn_i in function_list:
            self.edge_list.append((0, ftn_i))
            ftn_mark[ftn_i] = True
            if ftn_i not in self.nodes_w_edges:
                self.nodes_w_edges.append(ftn_i)

        self.nodes_wo_edges = list(set(self.nodes).difference(set(self.nodes_w_edges)))
        # ignore nodes that have no edges
        if len(self.nodes_wo_edges) > 0:
            for ftn in self.nodes_wo_edges:
                ftn_mark[ftn] = True

        # ---------------------------
        # at the depth 1, add edges from construrctor to functions in the function_list

        # set to the value of self.num_label
        for ftn_idx in function_list:
            self.fdg_3d_array[0, ftn_idx, 0] = self.num_label


        # ---------------------------
        # at the depth > 1
        for i in range(1, self.limit_maximum_depth):  # the maximum depth to 10
            # get function indices that can be reached from the immediate previous step
            ftn_reached = self.ftn_reached_at_a_depth(-1,False)

            array_i = -np.ones((1, self.num_ftn, self.num_ftn))
            # for each of these functions,get functions that depend on it.
            for idx_from in ftn_reached:
                write_sv_ary = self.sv_write[int(idx_from),:]
                # get sv indices that ftn (with index idx_from)  writes
                sv_indices_w = write_sv_ary[write_sv_ary >= 0]

                ftn_read=[]
                for sv_idx in sv_indices_w:
                    # get indices of functions that read label with index sv_idx
                    ftn_indices = self.sv_read[:, int(sv_idx)]
                    ftn_indices_r = np.where(ftn_indices >= 0)[0]
                    ftn_read+=list(ftn_indices_r) # save all functions reading sv_idx
                    if len(ftn_indices_r)>0:
                        for idx_to in ftn_indices_r:
                            if idx_to != idx_from:
                                array_i[0, int(idx_to), int(idx_from)] = sv_idx
                                if (int(idx_from),int(idx_to)) in edge_not_reached:
                                    edge_not_reached.remove((int(idx_from),int(idx_to)))

                # mark the functions reached
                for idx in ftn_read:
                    ftn_mark[idx]=True


            # ---------------------------
            # when no function is reached, stop
            comparison_0=array_i==-np.ones((1, self.num_ftn, self.num_ftn))
            if comparison_0.all():
                break
            comparison_1=array_i[0]==self.fdg_3d_array[-1]
            if comparison_1.all():
                break

            self.fdg_3d_array = np.concatenate((self.fdg_3d_array, array_i), axis=0)

            # when all functions are reached at least once,
            if all(ftn_mark):
                if self.depth_all_ftns_reached==0:
                    self.depth_all_ftns_reached = self.fdg_3d_array.shape[0]

            #when all edges are considered at least once, stop
            if len(edge_not_reached)==0:
                self.depth_all_edges_reached=self.fdg_3d_array.shape[0]
                break

        if self.depth_all_edges_reached==0:
            self.depth_all_ftns_reached=self.fdg_3d_array.shape[0]
            self.depth_all_edges_reached=self.fdg_3d_array.shape[0]
        # # assume that self.depth_all_ftns_reached<=self.depth_all_edges_reached
        # assert(self.depth_all_ftns_reached<=self.depth_all_edges_reached)












if __name__=='__main__':
    # if len(sys.argv)==5:
    #     file_path=sys.argv[1]+sys.argv[2]
    #     ftn_info= Function_info(file_path,sys.argv[3])
    #     functionsDict = ftn_info.functions_dict_slither()
    #     fdg=FDG(functionsDict )
    #     fdg.write_CSV(sys.argv[4])
    #     fdg.useful_data()
    # elif len(sys.argv)==3:
    #     write_csv_useful_data(sys.argv[1],sys.argv[2])
    # else:
    #     print("please provide three arguments: directory, solidity file, and contract name")
    #
    # # #==============
    # # file='/home/wei/PycharmProjects/Contracts/example/FDG/FDG_useful_data.txt'
    # # write_csv_useful_data(file, 'sss.csv')


    # #---------------------------------------------------------------------
    # # # get data needed to draw graph
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/example/contracts/ZetTokenMint.sol','ZetTokenMint')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/wqp_experiment_2021_april_2/contracts_exp/TrueCADMock.sol', 'Registry')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/wqp_experiment_2021_april_2/contracts_exp/AirDropContract.sol','AirDropContract')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/wqp_experiment_2021_april_2/contracts_exp/DocumentSigner.sol','DocumentSigner')
    # ftn_info = Function_info(
    #     '/home/wei/PycharmProjects/Contracts/wqp_experiment_2021_april_2/contracts_exp/play_me_quiz.sol',
    #     'play_me_quiz')

    ftn_info=Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')

    functionsDict=ftn_info.functions_dict_slither()


    # functionsDict={'f1': ['transferOwnership(address)', ['owner'], ['owner'], 'f2fde38b', 1], 'f2': ['transfer(address,uint256)', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb', 2], 'f3': ['transferFrom(address,address,uint256)', ['allowed', 'balances', 'mintingFinished'], ['allowed', 'balances'], '23b872dd', 3], 'f4': ['approve(address,uint256)', ['mintingFinished'], ['allowed'], '095ea7b3', 4], 'f5': ['increaseApproval(address,uint256)', [], ['allowed'], 'd73dd623', 5], 'f6': ['decreaseApproval(address,uint256)', [], ['allowed'], '66188463', 6], 'f7': ['setMinter(address)', ['owner'], ['minter'], 'fca3b5aa', 7], 'f8': ['mint(address,uint256)', ['balances', 'totalSupply', 'minter', 'mintingFinished'], ['balances', 'totalSupply'], '40c10f19', 8], 'f9': ['finishMinting()', ['minter'], ['mintingFinished'], '7d64bcb4', 9], 'f10': ['setDestroyer(address)', ['owner'], ['destroyer'], '6a7301b8', 10], 'f11': ['burn(uint256)', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68', 11], 'f0': ['constructor', [], ['owner'], '8afc3605', 0]}


    fdg=FDG(functionsDict)

    function_list=[1,5,7,6,10]

    fdg.build_fdg_3d_array(function_list)
    print('f')






    # # fdg = FDG(Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken').functions_dict_slither())
    # # print(fdg.matrix_fdg)
    # # print(fdg.depth_all_ftns_reached)
    #
    #
    # # #====================================
    # # functions_dict={'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'], 'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'], 'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'], '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'], 'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'], 'f6': ['decreaseApproval', [], ['allowed'], '66188463'], 'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'], 'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'], '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'], 'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'], 'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'], 'f0': ['constructor', [], ['owner'], '8afc3605']}
    # # fdg=FDG(functions_dict)
    # # print(fdg.depth_all_ftns_reached)
    #
    # #
    # # a=np.array([2,3,4,5])
    # # b=a[a>2]
    # # print(a)
    # # print(b)
    #
    # # print(set(fdg.ftn_reached_at_a_depth(1)).difference(set(fdg.ftn_reached_at_a_depth(2))))
    #
    # # fdg.write_CSV('/media/sf___share_vms/fdg_matrix.csv')








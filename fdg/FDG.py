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
        self.edges_no_f=[] # used for topological ordering of nodes in sequence generation
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
                    if ftn_from != ftn_to: # does not consider self data dependency
                        if ftn_from in self.graph.keys():
                            self.graph[ftn_from] += [ftn_to]
                        else:
                            self.graph[ftn_from] = [ftn_to]

                        if ftn_to not in self.nodes_w_DD_edges:
                            self.nodes_w_DD_edges.append(ftn_to)

                        # save edges
                        if 'f'+str(ftn_from)+",f"+str(ftn_to) not in self.edges.keys():
                            self.edges['f'+str(ftn_from)+",f"+str(ftn_to)]=self.index_to_label[sv_w_idx]
                        if [ftn_from, ftn_to] not in self.edges_no_f:
                            self.edges_no_f.append([ftn_from, ftn_to])

                # save nodes that have edges
                if ftn_from not in self.nodes_w_DD_edges:
                    self.nodes_w_DD_edges.append(ftn_from)

                # remove repeated element
                if ftn_from in self.graph.keys():
                    self.graph[ftn_from]=list(set(self.graph[ftn_from]))

        self.nodes_wo_DD_edges = list(set(self.nodes).difference(set(self.nodes_w_DD_edges)))
        self.nodes_wo_DD_edges=[item for item in self.nodes_wo_DD_edges if item not in [0,1]]

        self.depth_limit=5




















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














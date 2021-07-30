from copy import copy

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from fdg.funtion_info import Function_info
from fdg.sequence import Sequence

colors = ['orange', 'green', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'navy', 'blueviolet', 'purple',
          'magenta', 'crimson']


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
            self.index_to_ftn[value[4]]=value[0]
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

            if value[1]:
                for label in value[1]:
                    lbl_index=self.label_to_index[label]
                    read[lbl_index] = int(lbl_index)

            if value[2]:
                for label in value[2]:
                    lbl_index = self.label_to_index[label]
                    write[lbl_index] = int(lbl_index)

            if value[4]==0:
                self.sv_read[0,:]=-np.ones(self.num_label,dtype=int).T
                self.sv_write[0,:] = -np.ones(self.num_label,dtype=int).T
            else:
                self.sv_read[value[4],:]=read.T
                self.sv_write[value[4],:] = write.T

        # get all edges
        self.nodes = list(range(self.num_ftn))
        self.nodes_wo_edges=[]
        self.nodes_w_edges=[]
        self.start_nodes=[]
        self.graph = {}
        self.edge_list=[]

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
                        if ftn_to not in self.nodes_w_edges:
                            self.nodes_w_edges.append(ftn_to)
                        self.edge_list.append((ftn_from, ftn_to))
                # save nodes that have edges
                if ftn_from not in self.nodes_w_edges:
                    self.nodes_w_edges.append(ftn_from)

        # remove repeated element
        for ftn_from in self.graph.keys():
            self.graph[ftn_from]=list(set(self.graph[ftn_from]))


        self.fdg_3d_array=-np.ones((1, self.num_ftn, self.num_ftn))
        self.fdg_3d_array_new=-np.ones((1, self.num_ftn, self.num_ftn))

        # set the default value for depth_all_ftns_reached
        self.depth_all_ftns_reached = 0
        self.depth_all_edges_reached=0



    def build_fdg_3d_array(self,function_list:list):
        """
        stop build FDG when all edges are included
        :param function_list: [10,3,6]
        :return:
        """
        assert(len(function_list)>0)
        edge_not_reached = copy(self.edge_list)
        ftn_mark = [False] * self.num_ftn
        ftn_mark[0] = True  # for constructor,set it true

        self.graph[0]=function_list
        for ftn_i in function_list:
            self.edge_list.append((0,ftn_i))
            ftn_mark[ftn_i] = True
            if ftn_i not in self.nodes_w_edges:
                self.nodes_w_edges.append(ftn_i)

        self.nodes_wo_edges = list(set(self.nodes).difference(set(self.nodes_w_edges)))
        # ignore nodes that have no edges
        if len(self.nodes_wo_edges)>0:
            for ftn in self.nodes_wo_edges:
                ftn_mark[ftn]=True

        # ---------------------------
        # at the depth 1, add edges from construrctor to functions in the function_list

        # set to the value of self.num_label
        for ftn_idx in function_list:
            self.fdg_3d_array[0,ftn_idx,0]=self.num_label
            self.fdg_3d_array_new[0, ftn_idx, 0] = self.num_label
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
                    if len(ftn_indices_r)>0:
                        ftn_read += list(ftn_indices_r)  # save all functions reading sv_idx
                        for idx_to in ftn_indices_r:
                            if idx_to != idx_from:
                                array_i[0, int(idx_to), int(idx_from)] = sv_idx
                                # remove the edge reached
                                if (int(idx_from),int(idx_to)) in edge_not_reached:
                                    edge_not_reached.remove((int(idx_from),int(idx_to)))
                                # mark the node reached
                                if not ftn_mark[idx_to]:
                                    ftn_mark[idx_to]=True
            # stop conditions
            # ---------------------------
            # when no function is reached, stop
            comparison_0=array_i==-np.ones((1, self.num_ftn, self.num_ftn))
            if comparison_0.all():
                break
            # when the last lement is equal to the current element
            comparison_1=array_i[0]==self.fdg_3d_array[-1]
            if comparison_1.all():
                break

            self.fdg_3d_array = np.concatenate((self.fdg_3d_array, array_i), axis=0)

            # when all functions are reached at least once,save the depth value
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

        # assume that self.depth_all_ftns_reached<=self.depth_all_edges_reached
        # assert(self.depth_all_ftns_reached<=self.depth_all_edges_reached)


    def build_fdg_3d_array_new(self, function_pairs_dict: dict):
        """
        build new fdg from depth 2 to a depth all functions are reached
        :param function_pairs_dict:
        :return:
        """
        assert self.fdg_3d_array.shape[0] >= 1
        assert self.fdg_3d_array_new.shape[0] == 1
        for depth in range(1, self.fdg_3d_array.shape[0]):
            array_d = -np.ones((1, self.num_ftn, self.num_ftn))
            if depth in function_pairs_dict.keys():
                ftn_pairs = function_pairs_dict[depth]
                for (ftn_from, ftn_to) in ftn_pairs:
                    array_d[0, ftn_to, ftn_from] = self.fdg_3d_array[depth, ftn_to, ftn_from]
                self.fdg_3d_array_new = np.concatenate((self.fdg_3d_array_new, array_d), axis=0)
            else:
                array_d[0, :, :] = self.fdg_3d_array[depth, :, :]
            self.fdg_3d_array_new = np.concatenate((self.fdg_3d_array_new, array_d), axis=0)


    def ftn_reached_at_a_depth(self,depth,new_fdg_flag:bool)->list:
        """
        :param depth:
        :return: a list of function indices that reachable
        """

        if new_fdg_flag:
            ftn_reachable = self.fdg_3d_array_new[depth, :, :]
        else:
            ftn_reachable = self.fdg_3d_array[depth, :, :]

        row_indices=np.unique(np.nonzero(ftn_reachable>=0)[0])
        return row_indices

    def ftn_parents_reached_at_a_depth(self, depth) -> list:
        """
        :param depth:
        :return: a list of function indices that reachable
        """
        re = self.fdg_3d_array[depth, :, :]
        col_indices = np.unique(np.nonzero(re >= 0)[1])
        return col_indices

class fdg_model_1():
    def __init__(self, fdg:np.ndarray):
        self.fdg_3d_array=fdg
        self.graph_backward={}
        # graph_backward_parent_label={}
        for i in range(self.fdg_3d_array.shape[0] - 1, -1, -1):
            ftn_reached = self.ftn_reached_at_a_depth(i)
            for ftn_from in ftn_reached:
                ftn_parents_all = self.fdg_3d_array[i, ftn_from, :]
                ftn_parents = list(np.where(ftn_parents_all >= 0)[0])
                if ftn_from in ftn_parents:
                    ftn_parents.remove(ftn_from)
                if len(ftn_parents)==0:continue
                if ftn_from not in self.graph_backward.keys():
                    self.graph_backward[ftn_from] = list(ftn_parents)
                else:
                    self.graph_backward[ftn_from]+=list(ftn_parents)

        for key in self.graph_backward.keys():
            self.graph_backward[key]=list(set(self.graph_backward[key]))


    def ftn_reached_at_a_depth(self,depth)->list:
        """
        :param depth:
        :return: a list of function indices that reachable
        """
        ftn_reachable = self.fdg_3d_array[depth, :, :]
        row_indices=np.unique(np.nonzero(ftn_reachable>=0)[0])
        return row_indices

    def get_sequences(self, ftn_idx, graph: dict):
        def _get_sequences(graph: dict, function_idx: int, visited: list, path: list, seq: list):
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

        seq = []
        path = []
        visited = [False] * self.fdg_3d_array.shape[1]
        _get_sequences(graph, ftn_idx, visited, path, seq)
        return seq

    def draw_FDG_w_edge_label_backward(self, colors_list):
        nodes_list = ['f' + str(i) for i in range(self.fdg_3d_array.shape[1])]
        # get needed edge data from edges_label_dictt
        edges_list = []
        edges_label_dict = {}
        edges_label_list = []
        for i in range(self.fdg_3d_array.shape[0] - 1, -1, -1):
            ftn_reached = self.ftn_reached_at_a_depth(i)
            for ftn_to in ftn_reached:
                ftn_parents_all = self.fdg_3d_array[i, ftn_to, :]
                ftn_parents = np.where(ftn_parents_all >= 0)[0]
                if len(ftn_parents) == 0: continue
                for ftn_from in ftn_parents:
                    if ftn_from !=ftn_to:
                        e_tuple = ('f' + str(ftn_to), 'f' + str(ftn_from))
                        edges_list.append(e_tuple)
                        edges_label_dict[e_tuple] = self.fdg_3d_array[i,ftn_to,ftn_from]
                        edges_label_list.append(self.fdg_3d_array[i,ftn_to,ftn_from])

        # # ValueError: dictionary update sequence element #0 has length 1; 2 is required
        # # The problem is here tuple(e) a single item in the tuple, wants two items tuple(item1, item2)

        # create graph, add nodes and edges
        G = nx.MultiDiGraph()
        G.add_nodes_from(nodes_list)
        G.add_edges_from(edges_list)

        # about positioning nodes for the graph
        # pos = nx.layout.spectral_layout(G)
        # pos = nx.layout.planar_layout(G)
        pos = nx.layout.circular_layout(G)
        # pos = nx.layout.spring_layout(G)

        # get an color for each node
        nodes_color_list = []
        # only the start node had a different color
        for item in nodes_list:
            if item == 'f0':
                nodes_color_list.append('crimson')
                continue
            nodes_color_list.append('wheat')

        # assign a color to each unique edge lable
        edges_labels_unique = list(set(edges_label_list))
        edges_label_color_dict = {}

        for i in range(len(edges_labels_unique)):
            edges_label_color_dict[edges_labels_unique[i]] = colors_list[i]

        # for edges in the Graph, edges are ordered by item1,item2 in the edge tuple(item1,item2)
        # so, the edge colors should be based on the order of edges in G not in edges_list
        edges_color_list = []
        for edge in G.edges():
            label = edges_label_dict.get(edge)
            edges_color_list.append(edges_label_color_dict.get(label))

        # draw nodes
        nodes = nx.draw(G, pos, node_size=500, node_color=nodes_color_list, alpha=0.9,
                        labels={node: node for node in G.nodes()}, with_labels=True, font_color='black')

        # draw edges
        nx.draw_networkx_edges(G, pos, arrowstyle="->", arrowsize=20, edge_color=edges_color_list, width=2)

        # draw edge labels
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edges_label_dict, font_color='black')

        ax = plt.gca()
        ax.set_axis_off()
        plt.show()


    def draw_FDG_w_edge_label_forward(self, colors_list):
        nodes_list = ['f' + str(i) for i in range(self.fdg_3d_array.shape[1])]
        # get needed edge data from edges_label_dictt
        edges_list = []
        edges_label_dict = {}
        edges_label_list = []
        for i in range(self.fdg_3d_array.shape[0]):
            ftn_reached = self.ftn_reached_at_a_depth(i)
            for ftn_to in ftn_reached:
                ftn_parents_all = self.fdg_3d_array[i, ftn_to, :]
                ftn_parents = np.where(ftn_parents_all >= 0)[0]
                if len(ftn_parents) == 0: continue
                for ftn_from in ftn_parents:
                    if ftn_to != ftn_from:
                        e_tuple = ('f' + str(ftn_from), 'f' + str(ftn_to))
                        edges_list.append(e_tuple)
                        edges_label_dict[e_tuple] = self.fdg_3d_array[i, ftn_to, ftn_from]
                        edges_label_list.append(self.fdg_3d_array[i, ftn_to, ftn_from])



        # # ValueError: dictionary update sequence element #0 has length 1; 2 is required
        # # The problem is here tuple(e) a single item in the tuple, wants two items tuple(item1, item2)

        # create graph, add nodes and edges
        G = nx.MultiDiGraph()
        G.add_nodes_from(nodes_list)
        G.add_edges_from(edges_list)

        # about positioning nodes for the graph
        # pos = nx.layout.spectral_layout(G)
        # pos = nx.layout.planar_layout(G)
        pos = nx.layout.circular_layout(G)
        # pos = nx.layout.spring_layout(G)

        # get an color for each node
        nodes_color_list = []
        # only the start node had a different color
        for item in nodes_list:
            if item == 'f0':
                nodes_color_list.append('crimson')
                continue
            nodes_color_list.append('wheat')

        # assign a color to each unique edge lable
        edges_labels_unique = list(set(edges_label_list))
        edges_label_color_dict = {}

        for i in range(len(edges_labels_unique)):
            edges_label_color_dict[edges_labels_unique[i]] = colors_list[i]

        # for edges in the Graph, edges are ordered by item1,item2 in the edge tuple(item1,item2)
        # so, the edge colors should be based on the order of edges in G not in edges_list
        edges_color_list = []
        for edge in G.edges():
            label = edges_label_dict.get(edge)
            edges_color_list.append(edges_label_color_dict.get(label))

        # draw nodes
        nodes = nx.draw(G, pos, node_size=500, node_color=nodes_color_list, alpha=0.9,
                        labels={node: node for node in G.nodes()}, with_labels=True, font_color='black')

        # draw edges
        nx.draw_networkx_edges(G, pos, arrowstyle="->", arrowsize=20, edge_color=edges_color_list, width=2)

        # draw edge labels
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edges_label_dict, font_color='black')

        ax = plt.gca()
        ax.set_axis_off()
        plt.show()


if __name__ == '__main__':

    # # model Crowdsale.sol
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/Crowdsale.sol', 'Crowdsale')
    # functionInfo=ftn_info.functions_dict_slither()
    # num_nodes,edges,index_to_ftn=get_nodeNum_edges(functionInfo)
    # print(index_to_ftn)
    # fdg=fdg_model(num_nodes,edges,[1,2,3])
    # graph_b=fdg.generate_graph_backward()
    # ftn_idx=4

    # model HoloToken.sol
    ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    functionInfo = ftn_info.functions_dict_slither()

    fdg=FDG(functionInfo)
    fdg.build_fdg_3d_array([1, 5, 7, 6, 10])

    # fdg_model=fdg_model_1(fdg.fdg_3d_array)
    # # generate and print sequences for specified functions
    # ftn_idx=[11,2]
    # for ftn_i in ftn_idx:
    #     seq = fdg_model.get_sequences(ftn_i, fdg_model.graph_backward)
    #     seq.sort(key=len)
    #     print(f'sequence of {ftn_i}({len(seq)}):')
    #     for s in seq:
    #         print(f'\t{s}')
    # # fdg_model.draw_FDG_w_edge_label_backward(colors)
    # # fdg_model.draw_FDG_w_edge_label_forward(colors)


    edges_evaluated_successfully = {1:[(1, 10), (1, 7), (7, 9), (7, 8)],2:[(7, 9), (9, 4), (7, 8), (9, 3)],3:[(9,2)]}
    fdg.build_fdg_3d_array_new(edges_evaluated_successfully)
    # fdg_model=fdg_model_1(fdg.fdg_3d_array_new)
    # # generate and print sequences for specified functions
    # ftn_idx=[11,2]
    # for ftn_i in ftn_idx:
    #     seq = fdg_model.get_sequences(ftn_i, fdg_model.graph_backward)
    #     seq.sort(key=len)
    #     print(f'sequence of {ftn_i}({len(seq)}):')
    #     for s in seq:
    #         print(f'\t{s}')
    # fdg_model.draw_FDG_w_edge_label_backward(colors)
    # fdg_model.draw_FDG_w_edge_label_forward(colors)




    seq_object = Sequence(fdg, [11], 5)
    # get ready for sequence generation(parpare some data for sequence generation)
    seq_object.prepare_sequence_generation()

    sequences_dict = seq_object.get_sequences_by_level(2)
    print(sequences_dict)
    sequences_dict = seq_object.get_sequences_by_level(3)
    print(sequences_dict)

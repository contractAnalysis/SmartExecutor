# import numpy as np
# import matplotlib.pyplot as plt
# import networkx as nx
#
# from fdg.funtion_info import Function_info
#
# colors = ['orange','green','purple','brown','pink','gray','olive','cyan','navy','blueviolet','purple','magenta','crimson']
#
# def get_nodeNum_edges(functionInfo:dict):
#     ftn_to_index = {}
#     index_to_ftn = {}
#     labels=set()
#     for key, value in functionInfo.items():
#         # value: full_name,read_set, write_set,hash_value, index
#         ftn_to_index[value[0]] = value[4]
#         index_to_ftn[value[4]] = value[0]
#         labels =labels.union(value[1] + value[2])
#
#     labels = sorted(labels)
#     index_to_label={}
#     label_to_index={}
#     for idx, item in enumerate(labels):
#         index_to_label[idx] = item
#         label_to_index[item] = idx
#
#     num_ftn = len(list(functionInfo.keys()))
#     num_label = len(list(labels))
#
#     # create read, write matrices for all functions
#     sv_read = -np.ones((num_ftn, num_label), dtype=int)
#     sv_write = -np.ones((num_ftn, num_label), dtype=int)
#     for key, value in functionInfo.items():
#         read = -np.ones(num_label, dtype=int)
#         write = -np.ones(num_label, dtype=int)
#
#         if value[1]:
#             for label in value[1]:
#                 lbl_index = label_to_index[label]
#                 read[lbl_index] = int(lbl_index)
#
#         if value[2]:
#             for label in value[2]:
#                 lbl_index = label_to_index[label]
#                 write[lbl_index] = int(lbl_index)
#
#         if value[4] == 0:
#             sv_read[0, :] = -np.ones(num_label, dtype=int).T
#             sv_write[0, :] = -np.ones(num_label, dtype=int).T
#         else:
#             sv_read[value[4], :] = read.T
#             sv_write[value[4], :] = write.T
#
#     # get all edges
#     edge_list = []
#     # ignore constructor
#     for ftn_from in range(1, num_ftn):
#         sv_write__ = sv_write[ftn_from, :]
#         sv_write_idx = sv_write__[sv_write__ >= 0]  # indices of state variables written
#         if len(sv_write_idx) == 0:
#             continue
#         for sv_w_idx in sv_write_idx:
#             ftn_to = sv_read[:, sv_w_idx]
#             ftn_to_idx = np.where(ftn_to >= 0)[0]  # indices of functions reading sv_w_idx
#             if len(ftn_to_idx) == 0:
#                 continue
#             for ftn_to in ftn_to_idx:
#                 if ftn_from != ftn_to:
#                     edge_list.append((ftn_from, ftn_to,index_to_label[sv_w_idx]))
#     return num_ftn,edge_list,index_to_ftn
#
#
#
#
# class fdg_model():
#     def __init__(self, num_nodes:int, edge_list,start_nodes:list):
#         self.nodes=list(range(num_nodes))
#         self.edges=edge_list
#         self.start_nodes=start_nodes
#
#     def generate_graph_forward(self):
#         graph={}
#         graph[0]=self.start_nodes
#         for edge in self.edges:
#             p=edge[0]
#             c=edge[1]
#             if p not in graph.keys():
#                 graph[p]=[c]
#             else:
#                 graph[p]+=[c]
#         return graph
#
#     def generate_graph_backward(self):
#         graph={}
#         for ftn_idx in self.start_nodes:
#             graph[ftn_idx]=[0]
#         for edge in self.edges:
#             p=edge[0]
#             c=edge[1]
#             if c not in graph.keys():
#                 graph[c]=[p]
#             else:
#                 graph[c]+=[p]
#         return graph
#
#     def get_sequences(self,ftn_idx,graph:dict):
#         def _get_sequences(graph:dict, function_idx:int, visited:list, path:list, seq:list):
#             """
#             The general method to get seqence from function_idx to leaf nodes in the graph
#             :param graph: a dict
#             :param function_idx: the source
#             :param visited: [# of nodes]
#             :param path:[]
#             :param seq: save the resulted sequences
#
#             """
#             # checking all the visited nodes
#             visited[function_idx] = True
#             path.append(function_idx)
#
#             if function_idx not in graph:  # leaf nodes
#                 # se.append(list(np.copy(path))[1::][::-1])
#                 seq.append(list(np.copy(path)))
#             else:
#                 for i in graph[function_idx]:
#                     if visited[i] == False:
#                         _get_sequences(graph, i, visited, path, seq)
#             path.pop()
#             visited[function_idx] = False
#         seq=[]
#         path = []
#         visited = [False] * len(self.nodes)
#         _get_sequences(graph, ftn_idx, visited, path, seq)
#         return seq
#
#     def draw_FDG_w_edge_label(self, colors_list):
#         nodes_list=['f'+str(i) for i in self.nodes]
#         # get needed edge data from edges_label_dictt
#         edges_list = []
#         edges_label_dict = {}
#         edges_label_list = []
#
#         for edge in self.edges:
#             e_tuple=('f' + str(edge[0]),'f' + str(edge[1]))
#             edges_list.append(e_tuple)
#             edges_label_dict[e_tuple]=edge[2]
#             edges_label_list.append(edge[2])
#
#         for ftn_idx in self.start_nodes:
#             e_tuple=('f' + str(0),'f' + str(ftn_idx))
#             edges_list.append(e_tuple)
#             edges_label_dict[e_tuple]='none'
#             edges_label_list.append('none')
#
#         # # ValueError: dictionary update sequence element #0 has length 1; 2 is required
#         # # The problem is here tuple(e) a single item in the tuple, wants two items tuple(item1, item2)
#
#         # create graph, add nodes and edges
#         G = nx.MultiDiGraph()
#         G.add_nodes_from(nodes_list)
#         G.add_edges_from(edges_list)
#
#         # about positioning nodes for the graph
#         # pos = nx.layout.spectral_layout(G)
#         # pos = nx.layout.planar_layout(G)
#         pos = nx.layout.circular_layout(G)
#         # pos = nx.layout.spring_layout(G)
#
#         # get an color for each node
#         nodes_color_list = []
#         # only the start node had a different color
#         for item in nodes_list:
#             if item == 'f0':
#                 nodes_color_list.append('crimson')
#                 continue
#             nodes_color_list.append('wheat')
#
#         # assign a color to each unique edge lable
#         edges_labels_unique = list(set(edges_label_list))
#         edges_label_color_dict = {}
#
#         for i in range(len(edges_labels_unique)):
#             edges_label_color_dict[edges_labels_unique[i]] = colors_list[i]
#
#         # for edges in the Graph, edges are ordered by item1,item2 in the edge tuple(item1,item2)
#         # so, the edge colors should be based on the order of edges in G not in edges_list
#         edges_color_list = []
#         for edge in G.edges():
#             label = edges_label_dict.get(edge)
#             edges_color_list.append(edges_label_color_dict.get(label))
#
#         # draw nodes
#         nodes = nx.draw(G, pos, node_size=500, node_color=nodes_color_list, alpha=0.9,
#                         labels={node: node for node in G.nodes()}, with_labels=True, font_color='black')
#
#         # draw edges
#         nx.draw_networkx_edges(G, pos, arrowstyle="->", arrowsize=20, edge_color=edges_color_list, width=2)
#
#         # draw edge labels
#         nx.draw_networkx_edge_labels(G, pos, edge_labels=edges_label_dict, font_color='black')
#
#         ax = plt.gca()
#         ax.set_axis_off()
#         plt.show()
#
#
# if __name__=='__main__':
#
#     # # model Crowdsale.sol
#     # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/Crowdsale.sol', 'Crowdsale')
#     # functionInfo=ftn_info.functions_dict_slither()
#     # num_nodes,edges,index_to_ftn=get_nodeNum_edges(functionInfo)
#     # print(index_to_ftn)
#     # fdg=fdg_model(num_nodes,edges,[1,2,3])
#     # graph_b=fdg.generate_graph_backward()
#     # ftn_idx=4
#
#     # model HoloToken.sol
#     ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
#     functionInfo=ftn_info.functions_dict_slither()
#     num_nodes,edges,index_to_ftn=get_nodeNum_edges(functionInfo)
#     print(index_to_ftn)
#
#     fdg=fdg_model(num_nodes,edges,[1, 5, 7, 6, 10])
#     graph_b=fdg.generate_graph_backward()
#     ftn_idx=11
#     seq=fdg.get_sequences(ftn_idx,graph_b)
#     seq.sort(key=len)
#     print(f'sequence of {ftn_idx}({len(seq)}):')
#     for s in seq:
#         print(f'\t{s}')
#
#     ftn_idx = 2
#     seq = fdg.get_sequences(ftn_idx, graph_b)
#     seq.sort(key=len)
#     print(f'sequence of {ftn_idx}({len(seq)}):')
#     for s in seq:
#         print(f'\t{s}')
#     edges_kept=[(1,10),(1,7),(7,9),(7,8),(1,10),(1,10),(1,10),(1,10),(1,10),(1,10),(1,10)]
#
#     fdg.draw_FDG_w_edge_label(colors)
def merge_2_list(a, b) -> list:
    a1 = []  # temporarily save elements in a that are not matched with element in b
    b1 = []  # temporarily save elements in b that are not matched with element in a
    c = []  # save the combined result

    if len(b) == 0:
        c = a
        return c
    if len(a) == 0:
        c = b
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
                    else:  # both a and b still do not reach end, go to the next iteration
                        j += 1
                flag = False
                break
            else:  # do not match, save a[i] to a1
                a1.append(a[i])

        if flag:
            if j < len(b) - 1:
                a1 = []
                b1.append(b[j])
                j += 1
            else:
                c += a
                c.append(b[j])
                return c


def merge_fix_list_1(list_1, list_2):
    c = list_1
    j = 0
    while j < len(list_2):
        for i in range(len(list_1)):
            if j == len(list_2): break
            if list_1[i] == list_2[j]:
                j += 1

        if j < len(list_2):
            c += list_2[j:]
        break
    return c


def merge_fix_list_1_specified_lenth(list_1, list_2, specified_length):
    c = list_1[0:specified_length]
    j = 0
    while j < len(list_2):
        for i in range(specified_length):
            # if j == len(list_2): break
            # if i==len(list_2)-1
            # if list_1[i]==list_2[j]:
            #     j+=1
            if j == len(list_2):
                break;
            if list_1[i] != list_2[j]:
                break
            j += 1
        break
    a = list_1[specified_length - 1:]
    b = list_2[j:]
    re = merge_2_list(a, b)
    if re:
        c += re[1:]
    return c


if __name__ == '__main__':
    a = [1, 2, 5, 7]
    b = [4, 5]
    print(merge_2_list(a, b))

    # print(merge_fix_list_1([1,2,5,7],[4,5]))
    # print(merge_fix_list_1([1, 2, 5, 7], [2,8]))
    # print(merge_fix_list_1([1, 2, 5, 7], [2,7]))
    # print(merge_fix_list_1([1, 2, 5, 7], [2,1]))
    # print(merge_fix_list_1([1, 2, 5, 7], [1,3,7]))
    print('=====================')
    print(merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 5], 3))



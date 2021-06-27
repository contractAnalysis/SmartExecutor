import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from fdg.FDG_3d_array import *



def draw_FDG_w_edge_label(nodes_list,edges_dict,colors_list):

    # get needed edge data from edges_label_dictt
    edges_list = []
    edges_label_dict = {}
    edges_label_list = []

    for key, value in edges_dict.items():
        # edges_list.append(tuple(map(str,key.split(','))) )
        edge_tuple = tuple(key.split(','))
        edges_list.append(edge_tuple)
        edges_label_dict[edge_tuple] = value
        edges_label_list.append(value)
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

def draw_FDG_w_edge_label_node_label(functions_dict,nodes_list,edges_dict,colors_list):

    # get needed edge data from edges_label_dictt
    edges_list = []
    edges_label_dict = {}
    edges_label_list = []

    for key, value in edges_dict.items():
        # edges_list.append(tuple(map(str,key.split(','))) )
        edge_tuple = tuple(key.split(','))
        edges_list.append(edge_tuple)
        edges_label_dict[edge_tuple] = value
        edges_label_list.append(value)
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

    # get node name
    nodes_name_dict={}
    for node in G.nodes:
        nodes_name_dict[node]=functions_dict.get(node)[0]

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
    nodes = nx.draw(G, pos, node_size=5000, node_color=nodes_color_list, alpha=1,
                    labels={node: nodes_name_dict.get(node) for node in G.nodes()}, with_labels=True, font_color='black')

    # draw edges
    nx.draw_networkx_edges(G, pos, arrowstyle="->", arrowsize=20, edge_color=edges_color_list, width=2)

    # draw edge labels
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edges_label_dict, font_color='black')

    ax = plt.gca()
    ax.set_axis_off()
    plt.show()

def draw_FDG_wo_edge_label(nodes_list,edges_list):
    edges_list_1=[]
    for edge in edges_list:
        edge_tuple=tuple(str(edge).split(','))
        edges_list_1.append(edge_tuple)

    # create graph, add nodes and edges
    G = nx.MultiDiGraph()
    G.add_nodes_from(nodes_list)
    G.add_edges_from(edges_list_1)

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


    # draw nodes
    nodes = nx.draw(G, pos, node_size=500, node_color=nodes_color_list, alpha=0.9,
                    labels={node: node for node in G.nodes()}, with_labels=True, font_color='black',
                    arrowstyle="->", arrowsize=20,width=2
                    )

    ax = plt.gca()
    ax.set_axis_off()
    plt.show()

def draw_FDG_w_node_label(functions_dict,nodes_list,edges_dict):

    # get needed edge data from edges_label_dict
    edges_list = []
    for key, value in edges_dict.items():
        # edges_list.append(tuple(map(str,key.split(','))) )
        edge_tuple = tuple(key.split(','))
        edges_list.append(edge_tuple)

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

    # get node name
    nodes_name_dict={}
    for node in G.nodes:
        nodes_name_dict[node]=functions_dict.get(node)[0]

    # draw nodes
    nodes = nx.draw(G, pos, node_size=5000, node_color=nodes_color_list, alpha=0.9,
                    labels={node: nodes_name_dict.get(node) for node in G.nodes()}, with_labels=True, font_color='black',
                    arrowstyle="->", arrowsize=20,width=2
                    )

    ax = plt.gca()
    ax.set_axis_off()
    plt.show()


def get_nodes_edges_from_fdg(fdg:FDG):
    nodes=['f'+str(key) for key,_ in fdg.index_to_ftn.items()]
    edge_dict={}
    for i in range(fdg.fdg_3d_array.shape[0]):
        positions=np.nonzero(fdg.fdg_3d_array[i,:,:]>=0)
        for row,col in zip(positions[0],positions[1]):
            sv_index=fdg.fdg_3d_array[i,row,col]
            if sv_index==fdg.num_label:
                edge_dict['f'+str(col)+",f"+str(row)]='none'
            elif sv_index>=0:
                edge_dict['f'+str(col)+",f"+str(row)]=fdg.index_to_label[sv_index]
            else:
                edge_dict['f'+str(col)+",f"+str(row)]='does not exist'


    return nodes,edge_dict



# get edges
# for each function, check each of the state variables writen by it to see which functions read it.
# then, directed edges are created between the function and those functions.
def get_nodes_edges(dict):
    """
    # key: node pair, eg. f0,f1 (node name must be in format: f\d. eg., f0,f1,...,fn)
    # value: edge lable( state variable)the value of value of the dict must be an three-element array, and the last two elements are also arrays.
    :param dict:
        key: node name
        value: [ftn name, [ftn_list_write],[ftn_list_read]
    :return:
    """

    nodes_list=[]
    for key in dict.keys():
        nodes_list.append(key)

    edges_dict = {}
    # for each node ( or function)
    for i in range(len(dict)):
        # handle the case when the function does not read state variables
        sv_read=dict.get('f' + str(i))[1]

        if (len(sv_read)==0 and i!=0): # means no function writes a state variable that this function reads. ignore the case of constructor(f0)
            edges_dict['f0,f' + str(i)] = 'none'


        # handle ftn_list_read:
        sv_writen = dict.get('f' + str(i))[2]
        if len(sv_writen) > 0:
            for sv in sv_writen:
                for j in range(len(dict)):
                    if j != i:

                        if sv in dict.get('f' + str(j))[1]:
                            edges_dict['f' + str(i) + ',f' + str(j)] = sv
    return nodes_list,edges_dict




if __name__=='__main__':
    colors = ['orange','green','purple','brown','pink','gray','olive','cyan','navy','blueviolet','purple','magenta','crimson']

    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    # function_list = [1, 5, 7, 6, 10]
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/play_me_quiz.sol', 'play_me_quiz')
    # function_list = [1]
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/example_contracts/ZetTokenMint.sol', 'ZetTokenMint')
    # function_list = [1,2,4,6]
    ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/example_contracts/InvoiceCreator.sol', 'InvoiceCreator')
    function_list = [1,2]

    functionsDict = ftn_info.functions_dict_slither()
    # # visualize FDG based on function info
    # nodes,edges=get_nodes_edges(functionsDict)
    # print(f'number of nodes:{len(nodes)}')
    # print(f'number of edges:{len(edges)}')
    # draw_FDG_w_edge_label(nodes, edges,colors)

    #visualize FDG built from fdg_3d_array
    fdg=FDG(functionsDict)
    print(functionsDict)


    fdg.build_fdg_3d_array(function_list)
    print(f'{fdg.ftn_to_index}')

    nodes_,edges_=get_nodes_edges_from_fdg(fdg)
    print(f'number of nodes:{len(nodes_)}')
    print(f'nodes={nodes_}')
    print(f'number of edges:{len(edges_)}')
    print(f'edges={edges_}')

    draw_FDG_w_edge_label(nodes_, edges_,colors)



    # # visualize FDG built from fdg_3d_array, add
    # fdg=FDG(functionsDict)
    # nodes_1,edges_1=get_nodes_edges_from_fdg(fdg)
    # print(f'number of nodes:{len(nodes_1)}')
    # print(f'number of edges:{len(edges_1)}')
    # print(fdg.fdg_3d_array)
    # fdg_array1=fdg.fdg_3d_array
    # fdg.fdg_add(2,3)
    # nodes_2,edges_2=get_nodes_edges_from_fdg(fdg)
    # print(f'number of nodes:{len(nodes_2)}')
    # print(f'number of edges:{len(edges_2)}')
    # draw_FDG_w_edge_label(nodes_2, edges_2, colors)
    #
    #
    # # print(edges_1.keys())
    # # print(edges_2.keys())
    # # added_edges=list(set(edges_2).difference(set(edges_1)))
    # # print(f'added edges: {added_edges}')
    #
    # # fdg.fdg_delete(3)
    # # nodes_3,edges_3=get_nodes_edges_from_fdg(fdg)
    # # print(f'number of nodes:{len(nodes_3)}')
    # # print(f'number of edges:{len(edges_3)}')
    # # print(fdg.fdg_3d_array)
    # # fdg_array2=fdg.fdg_3d_array
    # # comparison=fdg_array1==fdg_array2
    # # if comparison.all():
    # #     print(f'True')
    # # else: print('False')
    # # # draw_FDG_w_edge_label(nodes_3, edges_3,colors)

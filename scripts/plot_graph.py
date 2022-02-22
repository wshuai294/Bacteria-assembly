import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import DBSCAN
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import re
from scipy.cluster.hierarchy import ward, fcluster
from scipy.spatial.distance import pdist
import pysam
import pandas as pd
import seaborn as sns
import networkx as nx
import sys

graph_file = sys.argv[1]
figure_file = sys.argv[2]
# graph_file = "/mnt/d/breakpoints/assembly/simulation/test/test.graph.txt"
# figure_file = "/mnt/d/breakpoints/assembly/simulation/test/test.graph.pdf"

def plot_net():
    g = open(graph_file)
    nodes_list = []
    edges_list = []
    first_partition_nodes = []
    first_partition_segs = {}
    for line in g:  
        array = line.strip().split()
        if array[0] == "SEG":
            continue
        seg1 = array[1]
        node1 = array[1] + array[2]
        seg2 = array[3]
        node2 = array[3] + array[4]
        if seg1 not in first_partition_segs and seg2 not in first_partition_segs:
            first_partition_segs[seg1] = 1
            first_partition_nodes.append(node1)
        elif seg1 in first_partition_segs:
            first_partition_nodes.append(node1)
        elif seg2 in first_partition_segs:
            first_partition_nodes.append(node2)           
        else:
            print ("both in the same side")
        nodes_list += [node1, node2]
        edges_list.append(tuple([node1, node2]))
    first_partition_nodes = sorted(first_partition_nodes)
    # DG = nx.Graph()
    DG = nx.DiGraph()
    DG.add_nodes_from(nodes_list)
    DG.nodes()
    print (len(edges_list), len(list(DG)), len(nodes_list))

    DG.add_edges_from(edges_list, weight = 0.1)
    
    pos = nx.spring_layout(DG, seed=8)        
    # nx.draw(DG, pos, with_labels=True, node_size=50, width=1, font_size=5)
    options = {
    'node_color': 'blue',
    'node_size': 50,
    'width': 1,
    'font_size': 5,
    'arrowstyle': '-|>',
    'arrowsize': 12,
}
    nx.draw_networkx(DG, pos, arrows=True, **options)
    # nx.draw(DG, 
    #     pos = nx.drawing.layout.bipartite_layout(DG, first_partition_nodes),  
    #     with_labels=True, 
    #     node_size=50, 
    #     width=1, 
    #     font_size=3)

    plt.savefig(figure_file)
    DG.clear()



def plot_net_BK():
    g = open(graph_file)
    nodes_list = []
    edges_list = []
    first_partition_nodes = []
    second_partition_nodes = []
    for line in g:  
        array = line.strip().split()
        node1 = array[0] + array[1]
        node2 = array[2] + array[3]
        if node1[-1] == "+":
            first_partition_nodes.append(node1)
        if node2[-1] == "+":
            first_partition_nodes.append(node2)
        nodes_list += [node1, node2]
        edges_list.append(tuple([node1, node2]))
    first_partition_nodes = sorted(first_partition_nodes)
    DG = nx.Graph()
    DG.add_nodes_from(nodes_list)
    DG.nodes()
    print (len(edges_list), len(list(DG)), len(nodes_list))

    DG.add_edges_from(edges_list, weight = 0.1)
    
    pos = nx.spring_layout(DG, seed=8)        
    # nx.draw(DG, pos, with_labels=True, node_size=50, width=1, font_size=5)
    nx.draw(DG, 
        pos = nx.drawing.layout.bipartite_layout(DG, first_partition_nodes),  
        with_labels=True, 
        node_size=50, 
        width=1, 
        font_size=5)

    plt.savefig(figure_file)
    DG.clear()

plot_net()
import pandas as pd
from itertools import combinations
import networkx as nx
import itertools
from collections import defaultdict
import time
import random
import seaborn as sns
import matplotlib.pyplot as plt

start_time = time.time()

def dictionary_index(orffile):

    orf_data=open(orffile, 'r')
    orfs=[]
    counter=0
    index={}
    for each_orf in orf_data:
        index[each_orf.strip()]=counter
        counter+=1

    inv_index = {v: k for k, v in index.iteritems()}
    return index, inv_index

dictionary_index("orf_yeast_list.csv")

orf_data=open("orf_yeast_list.csv", 'r')
orfs=[]
for each_orf in orf_data:
    orfs.append(each_orf.strip())

num_to_select=10
candidate_genes=random.sample(orfs, num_to_select)
#candidate_genes=random.seed(orfs, num_to_select)

def different_order_neighbours(biogrid_ppi, key_proteins, cutoff_level):
    
    index=dictionary_index("orf_yeast_list.csv")[0]
    inv_index=dictionary_index("orf_yeast_list.csv")[1]

###### PPI GRAPH #######
    file_ppi=open(biogrid_ppi, 'r')
    Graphtype=nx.Graph()
    G_ppi=nx.read_adjlist(file_ppi, create_using=Graphtype, nodetype=int)

    # MOTIF PATTERNS: TRIANGLES, SQUARES, PENTAGONS, HEXAGONS
    BFS_for_levels = nx.single_source_shortest_path_length(G_ppi, index[key_node], cutoff=cutoff_level)
    BFS_level1_neighbours=len([k for k,v in BFS_for_levels.items() if v == 1])
    BFS_level2_neighbours=len([k for k,v in BFS_for_levels.items() if v == 2])
    BFS_level3_neighbours=len([k for k,v in BFS_for_levels.items() if v == 3])
    BFS_level4_neighbours=len([k for k,v in BFS_for_levels.items() if v == 4])

    #print BFS_level1_neighbours, BFS_level2_neighbours, BFS_level3_neighbours, BFS_level4_neighbours

    return BFS_level1_neighbours, BFS_level2_neighbours, BFS_level3_neighbours, BFS_level4_neighbours

level1=[]
level2=[]
level3=[]
level4=[]
for key_node in candidate_genes:
    l1=different_order_neighbours("biogrid_all_ppi_new_graph.csv", key_node, 4)[0]
    l2=different_order_neighbours("biogrid_all_ppi_new_graph.csv", key_node, 4)[1]
    l3=different_order_neighbours("biogrid_all_ppi_new_graph.csv", key_node, 4)[2]
    l4=different_order_neighbours("biogrid_all_ppi_new_graph.csv", key_node, 4)[3]
    level1.append(l1)
    level2.append(l2)
    level3.append(l3)
    level4.append(l4)

avg_l1=sum(level1)/len(level1)
avg_l2=sum(level2)/len(level2)
avg_l3=sum(level3)/len(level3)
avg_l4=sum(level4)/len(level4)

sns.barplot(x=[1,2,3,4], y=[avg_l1, avg_l2, avg_l3, avg_l4], palette="Blues_d")
plt.xlabel("Level of neighbours (from the candidate node(level=0))")
plt.ylabel("Average number of nodes")
plt.title("Number of nodes in different layers of neighbours in PPI netowrk for a given protein")
#plt.xlim(xmin=0)
plt.show()
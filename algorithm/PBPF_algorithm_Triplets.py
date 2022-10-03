import pandas as pd
from itertools import combinations
import networkx as nx
import itertools
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import time
import numpy as np

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

def different_order_neighbours(biogrid_ppi, key_node, cutoff_level):
    
    index=dictionary_index("orf_yeast_list.csv")[0]
    inv_index=dictionary_index("orf_yeast_list.csv")[1]

###### PPI GRAPH #######
    file_ppi=open(biogrid_ppi, 'r')
    Graphtype=nx.Graph()
    G_ppi=nx.read_adjlist(file_ppi, create_using=Graphtype, nodetype=int)

    BFS_for_levels = nx.single_source_shortest_path_length(G_ppi, index[key_node.strip()], cutoff=cutoff_level)
    BFS_level1_neighbours=[k for k,v in BFS_for_levels.items() if v == 1]
    BFS_level2_neighbours=[k for k,v in BFS_for_levels.items() if v == 2]
    BFS_level3_neighbours=[k for k,v in BFS_for_levels.items() if v == 3]

    
    start_time3 = time.time()
    first_pairwise_combinations=list(itertools.permutations(BFS_level1_neighbours, 2))

    neighbours_connected_in_triangles = list(set(G_ppi.edges()).intersection(set(first_pairwise_combinations)))
    triangles_list=[[index[key_node.strip()]]+list(t) for t in neighbours_connected_in_triangles]
    number_triplets=len(triangles_list)
    end_time3=time.time()
    time_taken3=(end_time3 - start_time3)

    return number_triplets, time_taken3

triplets_number_time_outfile=open("output_triplets_transcription_factors.csv", "w")

n1=[]
t1=[]

for key_node in tumorigenesis_RP:
    level1_list_numbers=different_order_neighbours("biogrid_all_ppi_new_graph.csv", key_node, 3)[0]
    level1_list_time=different_order_neighbours("biogrid_all_ppi_new_graph.csv", key_node, 3)[1]
    #print squares_list_numbers, squares_list_time

    n1.append(level1_list_numbers)
    t1.append(level1_list_time)
    

data=pd.DataFrame({'number':n1,'Time-taken':t1})
#data=data[['Number of patterns', 'Time-taken']]
data.to_csv(triplets_number_time_outfile, index=False)


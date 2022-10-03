import pandas as pd
from itertools import combinations
import networkx as nx
import itertools
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import time

start_time3 = time.time()

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


def different_order_neighbours(biogrid_ppi, candidate_genes, cutoff_level):
    
    index=dictionary_index("orf_yeast_list.csv")[0]
    inv_index=dictionary_index("orf_yeast_list.csv")[1]

###### PPI GRAPH #######
    file_ppi=open(biogrid_ppi, 'r')
    Graphtype=nx.Graph()
    G_ppi=nx.read_adjlist(file_ppi, create_using=Graphtype, nodetype=int)

    module_list = pd.read_csv(candidate_genes, header=None)[0]
    key_node=module_list[1]

    BFS_for_levels = nx.single_source_shortest_path_length(G_ppi, index[key_node], cutoff=cutoff_level)
    BFS_level1_neighbours=[k for k,v in BFS_for_levels.items() if v == 1]
    BFS_level2_neighbours=[k for k,v in BFS_for_levels.items() if v == 2]
    BFS_level3_neighbours=[k for k,v in BFS_for_levels.items() if v == 3]

    keys = [0, 1, 2, 3]
    values = [[index[key_node]], BFS_level1_neighbours, BFS_level2_neighbours, BFS_level3_neighbours]
    dictionary = dict(zip(keys, values))

    first_pairwise_combinations=list(itertools.permutations(BFS_level1_neighbours, 2)) # for Triangles
    second_pairwise_combinations=list(itertools.permutations(BFS_level2_neighbours, 2)) # for Pentagons

    neighbours_connected_in_triangles = list(set(G_ppi.edges()).intersection(set(first_pairwise_combinations)))
    triangles_list=[[index[key_node]]+list(t) for t in neighbours_connected_in_triangles]

    end_time3=time.time()
    time_taken3=(end_time3 - start_time3)
    print("Triplets--- %s seconds ---" % (time_taken3))
    print("Number of Triplets: ", len(triangles_list))
    print("networkx triangles_list", nx.triangles(G_ppi, 275))  #782 (Benchmarked for triplets)

    start_time4 = time.time()

    squares_list=[]
    squares_no=[]
    for end_node in BFS_level2_neighbours:
        path_square=list(nx.all_simple_paths(G_ppi, source=index[key_node], target=end_node, cutoff=2))
        squares_no.append(1+2*(len(path_square)-2))
    number_of_squares=sum(squares_no)

    end_time4=time.time()
    time_taken4=(end_time4 - start_time4)
    print("Quadruplets--- %s seconds ---" % (time_taken4))
    print("Number of Quadruplets: ",number_of_squares)

    
    start_time5 = time.time()

    pentagons_nodes=[]
    end_edges_in_second_pairwise_combinations=set(G_ppi.edges()).intersection(second_pairwise_combinations)
    for end_node in end_edges_in_second_pairwise_combinations:
        paths_a = list(nx.all_simple_paths(G_ppi, source=index[key_node], target=end_node[0], cutoff=2))
        paths_b = list(nx.all_simple_paths(G_ppi, source=index[key_node], target=end_node[1], cutoff=2))
    print(paths_b)

    end_time5=time.time()
    time_taken5=(end_time5 - start_time5)
    print("Quintuplets--- %s seconds ---" % (time_taken5))
    print("Number of Quintuplets: ",len(pentagons_no))
    
    start_time6 = time.time()

    hexagons_list=[]
    hexagons_no=[]
    for end_node_level3 in BFS_level3_neighbours:
        path_hexagons=list(nx.all_simple_paths(G_ppi, source=index[key_node], target=end_node_level3, cutoff=3))
        hexagons_no.append(1+2*(len(path_hexagons)-2))
    number_of_hexagons=sum(hexagons_no)

    end_time6=time.time()
    time_taken6=(end_time6 - start_time6)
    print("Sextuplets--- %s seconds ---" % (time_taken6))
    print("Number of Sextuplets: ",number_of_hexagons)

    sns.barplot(x=[3,4,5,6], y=[time_taken3, time_taken4, time_taken5, time_taken6], palette="Blues_d")
    plt.xlabel("Size of connected-pattern")
    plt.ylabel("Time taken (in seconds)")
    plt.title("Time taken for searching different sized connected pattern in PPI-network")
    #plt.xlim(xmin=0)
    plt.show()
    

different_order_neighbours("biogrid_all_ppi_new_graph.csv", "chaperones.csv", 3)
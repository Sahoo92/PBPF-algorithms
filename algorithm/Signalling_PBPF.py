import pandas as pd
from itertools import combinations
import networkx as nx
import itertools
import time
from collections import Counter
from functools import partial
from operator import ne

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

def Interaction_Graph(biogrid_ppi, module):

	index=dictionary_index("orf_yeast_list.csv")[0]
	inv_index=dictionary_index("orf_yeast_list.csv")[1]

###### PPI GRAPH #######
	file_ppi=open(biogrid_ppi, 'r')
	Graphtype=nx.Graph()
	G_ppi=nx.read_adjlist(file_ppi, create_using=Graphtype, nodetype=int)
	print "# of nodes before: ", len(G_ppi.nodes())

	############ MODULE #################
	#print(index['YLR332C'])
	RP=open(module, "r")
	
	RP_list=[]
	for element in RP:
		RP_list.append(element)
	print RP_list
	RP_list.remove('YLR332C\n')

	module_list=[]
	for key_node in RP_list:
		#print key_node
		degree_key=G_ppi.degree(index[key_node.strip()])
		#print degree_key
		if(degree_key<400):
			module_list.append(index[key_node.strip()])
	#print module_list


	return G_ppi, module_list

Interaction_Graph("biogrid_all_ppi_new_graph.csv", "Yeast_signal_transduction.csv")

def motifs_in_modules(biogrid_ppi, key_node):

	index=dictionary_index("orf_yeast_list.csv")[0]
	inv_index=dictionary_index("orf_yeast_list.csv")[1]
	G_ppi=Interaction_Graph("biogrid_all_ppi_new_graph.csv", "Yeast_signal_transduction.csv")[0]
	
	highly_connected_nodes=[]
	for node_i in G_ppi.nodes():
		degree_i=G_ppi.degree(node_i)
		if(degree_i>400):
			highly_connected_nodes.append(node_i)
	#print "5555: ", 5555 in highly_connected_nodes

	G_ppi.remove_nodes_from(highly_connected_nodes)
	print "# of nodes after: ", len(G_ppi.nodes())

	print "# of nodes deleted: ", len(highly_connected_nodes) 
##### TRIANGLES ########

	start_time3 = time.time()
	#print "5555: ", 5555 in highly_connected_nodes
	print "*******************:", G_ppi.nodes()
	first_neighbours=list(G_ppi.neighbors(key_node))
	first_pairwise_combinations=list(itertools.permutations(first_neighbours, 2))
	
	neighbours_connected_in_triangles = list(set(G_ppi.edges()).intersection(set(first_pairwise_combinations)))
	triangles_list=[[key_node]+list(t) for t in neighbours_connected_in_triangles]

	end_time3=time.time()
	time_taken3=(end_time3 - start_time3)

	print inv_index[key_node], len(triangles_list)

##### UNIQUE LEVEL-1 PAIRWISE NODES  ################
	sorted_first_pairwise_combinations=[sorted(list(elem)) for elem in first_pairwise_combinations] 
	dict_first_pairwise=Counter([tuple(i) for i in sorted_first_pairwise_combinations])
	Output_dict_first_pairwise=pd.DataFrame(data={'list':list(dict_first_pairwise.keys()), 'count': list(dict_first_pairwise.values())})
	unique_first_pairwise_combinations_list_updated=Output_dict_first_pairwise['list']
    
##### SQUARES #########
	start_time4 = time.time()
	squares_list=[]

	for elements_squares in unique_first_pairwise_combinations_list_updated:
		n1=[n for n in G_ppi.neighbors(elements_squares[0])]
		n2=[n for n in G_ppi.neighbors(elements_squares[1])]
		n1.remove(key_node)
		n2.remove(key_node)	
		if(set(n1).intersection(set(n2))):
			for elements_square_list in set(n1).intersection(set(n2)):
				squares=[key_node, elements_squares[0], elements_square_list, elements_squares[1]]
				squares_list.append(sorted(squares))

	end_time4=time.time()
	time_taken4=(end_time4 - start_time4)

	print inv_index[key_node], len(squares_list)
              
##### PENTAGONS #########
	pentagons_list=[]
	second_neighbours1=[]
	second_neighbours2=[]
	second_pairwise_combinations={}
	for elements_pentagons in unique_first_pairwise_combinations_list_updated:
		n11=[n for n in G_ppi.neighbors(elements_pentagons[0])]
		n22=[n for n in G_ppi.neighbors(elements_pentagons[1])]
		second_neighbours1.append(n11)
		second_neighbours2.append(n22)
		n11.remove(key_node)
		n22.remove(key_node)

		second_pairwise1=list(itertools.product(n11,n22))
		second_pairwise2=list(itertools.product(n22,n11))
		second_pairwise_combinations[elements_pentagons]=second_pairwise1+second_pairwise2

		start_time5 = time.time()

		third_pairwise_combinations={}    
		third_pairwise_combinations_for_hexagons={}

	for keys in second_pairwise_combinations:
		values_key=sorted(second_pairwise_combinations[keys])
		list_values_keys=[list(elem) for elem in values_key]

		result1=list(filter(lambda x:x[0]==keys[0], values_key))
		result2=list(filter(lambda x:x[0]==keys[1], values_key))
		result3=list(filter(lambda x:x[1]==keys[0], values_key))
		result4=list(filter(lambda x:x[1]==keys[1], values_key))
		result5=list(filter(lambda x:x[0]==x[1], values_key))

		list_result1=[list(elem1) for elem1 in result1]
		list_result2=[list(elem2) for elem2 in result2]
		list_result3=[list(elem3) for elem3 in result3]
		list_result4=[list(elem4) for elem4 in result4]
		list_result5=[list(elem5) for elem5 in result5]
		all_list=list_result1+list_result2+list_result3+list_result4+list_result5
        
        tuple_all_list=[tuple(y) for y in all_list]

        l2_combinations=set(values_key)-set(tuple_all_list)
        third_pairwise_combinations[keys]=l2_combinations
        
        ################################ FOR HEXAGON #############################################
        unique_l2_combinations=set()
        test_unique_l2_combinations=[unique_l2_combinations.add((a,b)) for (a,b) in list(l2_combinations) if (a,b) and (b,a) not in unique_l2_combinations]
        third_pairwise_combinations_for_hexagons[keys]=list(unique_l2_combinations)
        ##########################################################################################
        
        neighbours_connected_in_pentagons = list(set(G_ppi.edges()).intersection(l2_combinations))
        list_neighbours_connected_in_pentagons=[list(a) for a in neighbours_connected_in_pentagons]
        if(len(list_neighbours_connected_in_pentagons)>0):
            for ii in list_neighbours_connected_in_pentagons:
                pentagons=[key_node, list(keys), ii]
                pentagons_list.append(pentagons)


        end_time5=time.time()
    	time_taken5=(end_time5 - start_time5)

    	print inv_index[key_node], len(pentagons_list)

####### HEXAGONS ##########
	start_time6 = time.time()

	hexagons_list=[]
	for keys_l3, v3 in third_pairwise_combinations_for_hexagons.iteritems():
		list_keysl3=list(keys_l3)
		for vv3 in v3:
			n_elements1=[n for n in G_ppi.neighbors(vv3[0])]
        	n_elements2=[n for n in G_ppi.neighbors(vv3[1])]

        	n_elements1_updated=list(filter(lambda num: num == 0, n_elements1))
        	n_elements2_updated=list(filter(lambda num: num == 0, n_elements2))
    		n_elements11=list(filter(lambda num: num == keys_l3[0], n_elements1))
    		n_elements22=list(filter(lambda num: num == keys_l3[0], n_elements2))
    		n_elements111=list(filter(lambda num: num == keys_l3[1], n_elements1))
    		n_elements222=list(filter(lambda num: num == keys_l3[1], n_elements2))

    		a1_list_l3=n_elements1_updated+n_elements11+n_elements111
    		b1_list_l3=n_elements2_updated+n_elements22+n_elements222

    		tuple_a1_list_l3=tuple(a1_list_l3)
    		tuple_b1_list_l3=tuple(b1_list_l3)

    		a1_l3_combinations=set(n_elements1)-set(tuple_a1_list_l3)
    		b1_l3_combinations=set(n_elements2)-set(tuple_b1_list_l3)

    		ss=a1_l3_combinations.intersection(b1_l3_combinations)

    		if(len(ss)>0):
    			for iii in ss:
    				hexagons=[key_node, keys_l3, vv3, iii]
    				hexagons_list.append(hexagons)
	end_time6=time.time()
	time_taken6=(end_time6 - start_time6)
	print inv_index[key_node], len(hexagons_list)

	return triangles_list, time_taken3, squares_list, time_taken4, pentagons_list, time_taken5, hexagons_list, time_taken6

inv_index=dictionary_index("orf_yeast_list.csv")[1]

Ribogenesis_PBPF_Motifs_Numbers_Time=open("Signalling_PBPF_Motifs_Numbers_Time.csv", "w")
module_list=Interaction_Graph("biogrid_all_ppi_new_graph.csv", "Yeast_signal_transduction.csv")[1]


n3=[]
t3=[]
n4=[]
t4=[]
n5=[]
t5=[]
n6=[]
t6=[]
protein=[]
for key_node in module_list:
	print key_node
	triangles_list=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[0]
	triangles_time=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[1]
	squares_list=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[2]
	squares_time=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[3]
	pentagons_list=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[4]
	pentagons_time=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[5]
	hexagons_list=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[6]
	hexagons_time=motifs_in_modules("biogrid_all_ppi_new_graph.csv", key_node)[7]

	print("Triplets: ", triangles_list, len(triangles_list), triangles_time)
	print("Quadruplets: ", squares_list, len(squares_list), squares_time)
	print("Quintuplets: ", pentagons_list, len(pentagons_list), pentagons_time)
	print("Sextuplets: ", hexagons_list, len(hexagons_list), hexagons_time)

	protein.append(inv_index[key_node])
	n3.append(len(triangles_list))
	t3.append(triangles_time)
	n4.append(len(squares_list))
	t4.append(squares_time)
	n5.append(len(pentagons_list))
	t5.append(pentagons_time)
	n6.append(len(hexagons_list))
	t6.append(hexagons_time)


data=pd.DataFrame({'protein': protein,'n3':n3,'t3':t3,'n4':n4,'t4':t4,'n5':n5,'t5':t5,'n6':n6,'t6':t6})
#data=data[['Number of patterns', 'Time-taken']]
data.to_csv(Signalling_PBPF_Motifs_Numbers_Time, index=False)


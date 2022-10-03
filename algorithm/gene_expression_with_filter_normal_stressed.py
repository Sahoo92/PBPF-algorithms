import pandas as pd
from itertools import combinations
import networkx as nx
import itertools
import time
from collections import Counter
from functools import partial
from operator import ne
import json

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

def expression_data_coregulation(yeast_expression_data):


	index=dictionary_index("orf_yeast_list.csv")[0]
	inv_index=dictionary_index("orf_yeast_list.csv")[1]


	#################### Modification to keep only the list of genes present in the indexing list #################

	list_of_genes=[]
	for g_index, genes in inv_index.iteritems():
		list_of_genes.append(genes)

	################################################################################################################

	data=pd.read_csv("complete_dataset.txt", sep="\t")
	normal=data["steady state 25 dec C ct-2"]
	normal_expression_values=normal.values.tolist()
	
	stressed=data["heat shock 25 to 37, 20 minutes"]
	stressed_expression_values=stressed.values.tolist()
	

	### NORMAL CONDITION ########
	gene_list_normal=[]

	for n1 in normal:
		#genes_nn=data.loc[normal == n1, 'UID']
		genes_nn1=data.loc[(normal == n1) & (data['steady state 25 dec C ct-2'] >2), 'UID']
		genes_nn2=data.loc[(normal == n1) & (data['steady state 25 dec C ct-2'] <-2), 'UID']
		frames1=[genes_nn1, genes_nn2]
		genes_nn=pd.concat(frames1)
		gene_list_normal.append(list(genes_nn))
	
	modified_normal_list = [item for sublist in gene_list_normal for item in sublist]

	normal_indexed_list=[]
	for i in set(modified_normal_list):
		if(i in list_of_genes):
			ind=index[i]
			normal_indexed_list.append(ind)

	
	### STRESSED CONDITION ########
	
	gene_list_stressed=[]

	for p2 in stressed:
		#genes_ps=data.loc[stressed == p2, 'UID']
		genes_ps1=data.loc[(stressed == p2) & (data['heat shock 25 to 37, 20 minutes'] >2), 'UID']
		genes_ps2=data.loc[(stressed == p2) & (data['heat shock 25 to 37, 20 minutes'] <-2), 'UID']
		frames2=[genes_ps1, genes_ps2]
		genes_ps=pd.concat(frames2)
		gene_list_stressed.append(list(genes_ps))

	modified_stressed_list = [item for sublist in gene_list_stressed for item in sublist]
	
	stressed_indexed_list=[]
	for j in set(modified_stressed_list):
		if(j in list_of_genes):
			indi=index[j]
			stressed_indexed_list.append(indi)

	print "normal", len(normal_indexed_list), "stressed", len(stressed_indexed_list)

	return normal_indexed_list, stressed_indexed_list, normal_expression_values, stressed_expression_values
 
expression_data_coregulation("complete_dataset.txt")


def gene_expression_integrated_motifs(motif):

	data=pd.read_csv("complete_dataset.txt", sep="\t")
	normal=data["steady state 25 dec C ct-2"]
	stressed=data["heat shock 25 to 37, 20 minutes"]

	#normal_gene_list=expression_data_coregulation("complete_dataset.txt")[0]
	stressed_gene_list=expression_data_coregulation("complete_dataset.txt")[1]
	#normal_expression_values=expression_data_coregulation("complete_dataset.txt")[2] 
	#stressed_expression_values=expression_data_coregulation("complete_dataset.txt")[3]

	motif_data=pd.read_csv(motif)
	
	############## TRIPLETS #######################################################
	Total_triplets={}
	for protein_names1 in motif_data["protein"]:
		triplets=motif_data.loc[motif_data["protein"] == protein_names1, 'n3']
		for trps in triplets:
			Total_triplets[protein_names1]=trps

	#print type(Total_triplets)

	############# SQUARES #########################################################
	Total_quadruplets={}
	for protein_names2 in motif_data["protein"]:
		quadruplets=motif_data.loc[motif_data["protein"] == protein_names2, 'n4']
		for quads in quadruplets:
			Total_quadruplets[protein_names2]=quads

	############## PENTAGONS ######################################################
	Total_quantuplets={}
	for protein_names3 in motif_data["protein"]:
		quantuplets=motif_data.loc[motif_data["protein"] == protein_names3, 'n5']
		for quants in quantuplets:
			Total_quantuplets[protein_names3]=quants

	############## HEXAGONS #######################################################
	Total_sextuplets={}
	for protein_names4 in motif_data["protein"]:
		sextuplets=motif_data.loc[motif_data["protein"] == protein_names4, 'n6']
		for hexs in sextuplets:
			Total_sextuplets[protein_names4]=hexs
	################################################################################	
	
	PI_triplets_total=[]
	PI_GI_triplets_total=[]
	PI_GE_Normal_triplets_total=[]
	PI_GE_Stressed_triplets_total=[]
	PI_GE_Stressed_triplets={}
	PI_GE_stressed_values=[]

	for k1 in Total_triplets:
		v1=Total_triplets[k1]
		print v1, type(v1)
		PI_triplets_total.append(len(v1))
		vv1=json.loads(v1)
		for vvv1 in vv1:
			#print("vvv1[]", vvv1[0])
			result1 = all(elem in stressed_gene_list for elem in vvv1)
			#s1=data.loc[vv1 == vvv1[0], 'steady state 25 dec C ct-2']
			#s2=data.loc[vv1 == vvv1[1], 'steady state 25 dec C ct-2']
			#s3=data.loc[vv1 == vvv1[2], 'steady state 25 dec C ct-2']
			#s123=sum(s1,s2,s3)
			if(result1==True): 
				PI_GE_Stressed_triplets_total.append(vvv1)
				PI_GE_Stressed_triplets[vvv1[0]]=vv1
				#PI_GE_stressed_values.append(s123)
	#print("PI_GE_Stressed_triplets: ", PI_GE_Stressed_triplets)
		#print("PI_GE_Stressed_triplets of ", vvv1[0], len(PI_GE_Stressed_triplets))
	#print(len(PI_GE_Stressed_triplets))
	#print ("PI_GE_Stressed_triplets_total: ", len(PI_GE_Stressed_triplets_total))
	#print ("PI_GE_stressed_values: ", PI_GE_stressed_values)

	for kk1 in PI_GE_Stressed_triplets:
		vkk1=PI_GE_Stressed_triplets[kk1]
		print len(vkk1)

	PI_quadruplets_total=[]
	PI_GI_quadruplets_total=[]
	PI_GE_Normal_quadruplets_total=[]
	PI_GE_Stressed_quadruplets_total=[]
	PI_GE_stressed_quadruplets={}

	for k2 in Total_quadruplets:
		v2=Total_quadruplets[k2]
		PI_quadruplets_total.append(len(v2))
		vv2=json.loads(v2)
		for vvv2 in vv2:
			result2 = all(elem in stressed_gene_list for elem in vvv2)
			if(result2==True): 
				PI_GE_Stressed_quadruplets_total.append(vvv2)
				PI_GE_stressed_quadruplets[vvv2[0]]=vv2

		#print ("PI_GE_stressed_quadruplets of ",vvv2[0], len(PI_GE_stressed_quadruplets))
	#print sum(PI_squares_total)
	#print("PI_GE_Stressed_quadruplets_total :", len(PI_GE_Stressed_quadruplets_total))

	for kk2 in PI_GE_stressed_quadruplets:
		vkk2=PI_GE_stressed_quadruplets[kk2]
		print "square", len(vkk2)

	PI_quantuplets_total=[]
	PI_GI_quantuplets_total=[]
	PI_GE_Normal_quantuplets_total=[]
	PI_GE_Stressed_quantuplets_total=[]
	PI_GE_stressed_quantuplets={}

	for k3 in Total_quantuplets:
		v3=Total_quantuplets[k3]
		PI_quantuplets_total.append(len(v3))
		vv3=json.loads(v3)
		for vvv3 in vv3:
			result3 = all(elem in stressed_gene_list for elem in vvv3)
			#PI_GE_stressed_quantuplets[vvv3[0]]=0
			if(result3==True): 
				PI_GE_Stressed_quantuplets_total.append(vvv3)
				PI_GE_stressed_quantuplets[vvv3[0]]=vv3

		#print ("PI_GE_stressed_quantuplets of ",vvv3[0], len(PI_GE_stressed_quantuplets))
	#print("PI_GE_Stressed_quantuplets_total: ", len(PI_GE_Stressed_quantuplets_total))
	#print sum(PI_pentagons_total)
	for kk3 in PI_GE_stressed_quantuplets:
		vkk3=PI_GE_stressed_quantuplets[kk3]
		print "PENTAGONS", len(vkk3)


	PI_hexagons_total=[]
	PI_GI_hexagons_total=[]
	PI_GE_Normal_hexagons_total=[]
	PI_GE_Stressed_hexagons_total=[]
	PI_GE_stressed_hexagons={}

	for k4 in Total_sextuplets:
		v4=Total_sextuplets[k4]
		PI_hexagons_total.append(len(v4))
		vv4=json.loads(v4)
		for vvv4 in vv4:
			result4 = all(elem in stressed_gene_list for elem in vvv4)
			if(result4==True): 
				PI_GE_Stressed_hexagons_total.append(vvv4)
				PI_GE_stressed_hexagons[vvv4[0]]=vv4

	for kk4 in PI_GE_stressed_hexagons:
		vkk4=PI_GE_stressed_hexagons[kk4]
		print "HEXAGONS", len(vkk4)

		#print ("PI_GE_stressed_hexagons of ", vvv4[0], len(PI_GE_stressed_hexagons))
	#print sum(PI_hexagons_total)
	#print("PI_GE_Stressed_hexagons_total :", len(PI_GE_Stressed_hexagons_total))

#gene_expression_integrated_motifs("PI_GI_Signalling_triplets.csv", "PI_GI_Signalling_quadruplets.csv", "PI_GI_Signalling_quintuplets.csv")
gene_expression_integrated_motifs("Signalling_degree_below300_PBPF_Motifs.csv")



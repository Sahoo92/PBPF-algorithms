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

    print inv_index[5770]
    
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
		genes_nn1=data.loc[(normal == n1) & (data['steady state 25 dec C ct-2'] >1), 'UID']
		genes_nn2=data.loc[(normal == n1) & (data['steady state 25 dec C ct-2'] <-1), 'UID']
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
		genes_ps1=data.loc[(stressed == p2) & (data['heat shock 25 to 37, 20 minutes'] >0.5), 'UID']
		genes_ps2=data.loc[(stressed == p2) & (data['heat shock 25 to 37, 20 minutes'] <-0.5), 'UID']
		frames2=[genes_ps1, genes_ps2]
		genes_ps=pd.concat(frames2)
		gene_list_stressed.append(list(genes_ps))

	modified_stressed_list = [item for sublist in gene_list_stressed for item in sublist]
	
	stressed_indexed_list=[]
	for j in set(modified_stressed_list):
		if(j in list_of_genes):
			indi=index[j]
			stressed_indexed_list.append(indi)

	print inv_index[5770]

	#print "normal", len(normal_indexed_list), "stressed", len(stressed_indexed_list)

	return normal_indexed_list, stressed_indexed_list, normal_expression_values, stressed_expression_values
 
expression_data_coregulation("complete_dataset.txt")


def gene_expression_integrated_motifs(motif1, motif2, motif3):

	inv_index=dictionary_index("orf_yeast_list.csv")[1]


	data=pd.read_csv("complete_dataset.txt", sep="\t")
	normal=data["steady state 25 dec C ct-2"]
	stressed=data["heat shock 25 to 37, 20 minutes"]

	normal_gene_list=expression_data_coregulation("complete_dataset.txt")[0]
	stressed_gene_list=expression_data_coregulation("complete_dataset.txt")[1]
	#normal_expression_values=expression_data_coregulation("complete_dataset.txt")[2] 
	#stressed_expression_values=expression_data_coregulation("complete_dataset.txt")[3]

	motif_data1=pd.read_csv(motif1)
	motif_data2=pd.read_csv(motif2)
	motif_data3=pd.read_csv(motif3)
	
	print inv_index[5770]
	############## TRIPLETS #######################################################
	Total_triplets={}
	for protein_names1 in motif_data1["protein"]:
		triplets=motif_data1.loc[motif_data1["protein"] == protein_names1, 'n3']
		for trps in triplets:
			Total_triplets[protein_names1]=trps

	#print (Total_triplets)

	############# SQUARES #########################################################
	Total_quadruplets={}
	for protein_names2 in motif_data2["protein"]:
		quadruplets=motif_data2.loc[motif_data2["protein"] == protein_names2, 'n3']
		for quads in quadruplets:
			Total_quadruplets[protein_names2]=quads

	############## PENTAGONS ######################################################
	Total_quantuplets={}
	for protein_names3 in motif_data3["protein"]:
		quantuplets=motif_data3.loc[motif_data3["protein"] == protein_names3, 'n3']
		for quants in quantuplets:
			Total_quantuplets[protein_names3]=quants

	
	################################################################################	
	
	PI_triplets_total=[]
	PI_GI_triplets_total=[]
	PI_GE_Normal_triplets_total=[]
	PI_GE_Stressed_triplets_total=[]
	PI_GE_Stressed_triplets={}
	PI_GE_stressed_values=[]

	stressed_list_proteins=[]

	for k1 in Total_triplets:
		v1=Total_triplets[k1]
		if(len(v1)>1):
			vv1=json.loads(v1)		
			for vvv1 in vv1:
				result1 = all(elem in stressed_gene_list for elem in vvv1)
				if(result1==True): 
					PI_GE_Stressed_triplets_total.append(vvv1)
					PI_GE_Stressed_triplets[vvv1[0]]=vv1
				
	#print PI_GE_Stressed_triplets

	for kk1 in PI_GE_Stressed_triplets:
		vkk1=PI_GE_Stressed_triplets[kk1]
		flat_list_vkk1 = [item for sublist in vkk1 for item in sublist]
		for vkk1_elem in flat_list_vkk1:
			stressed_list_proteins.append(inv_index[vkk1_elem])

			
	vkk1_elem_list=set(stressed_list_proteins)
	#for i in vkk1_elem_list:
		#print i

	PI_quadruplets_total=[]
	PI_GI_quadruplets_total=[]
	PI_GE_Normal_quadruplets_total=[]
	PI_GE_Stressed_quadruplets_total=[]
	PI_GE_stressed_quadruplets={}

	for k2 in Total_quadruplets:
		v2=Total_quadruplets[k2]
		if(len(v2)>1):
			vv2=json.loads(v2)
			for vvv2 in v2:
				result2 = all(elem in stressed_gene_list for elem in vvv2)
				if(result2==True): 
					PI_GE_Stressed_quadruplets_total.append(vvv2)
					PI_GE_stressed_quadruplets[vvv2[0]]=len(vv2)

	#print PI_GE_stressed_quadruplets


	for kk2 in PI_GE_stressed_quadruplets:
		vkk2=PI_GE_stressed_quadruplets[kk2]
		#print vkk2
	stressed_list_proteins_pent=[]
	PI_quantuplets_total=[]
	PI_GI_quantuplets_total=[]
	PI_GE_Normal_quantuplets_total=[]
	PI_GE_Stressed_quantuplets_total=[]
	PI_GE_stressed_quantuplets={}

	for k3 in Total_quantuplets:
		v3=Total_quantuplets[k3]
		if(len(v3)>1):
			vv3=json.loads(v3)
			for vvv3 in vv3:
				result3 = all(elem in stressed_gene_list for elem in vvv3)
			
				if(result3==True): 
					PI_GE_Stressed_quantuplets_total.append(vvv3)
					PI_GE_stressed_quantuplets[vvv3[0]]=vv3

	#print PI_GE_stressed_quantuplets


	for kk3 in PI_GE_stressed_quantuplets:
		vkk3=PI_GE_stressed_quantuplets[kk3]
		flat_list_vkk3 = [item for sublist in vkk3 for item in sublist]
		for vkk3_elem in flat_list_vkk3:
			stressed_list_proteins_pent.append(inv_index[vkk3_elem])
			
	vkk3_elem_list=set(stressed_list_proteins_pent)
	#for j in vkk3_elem_list:
	#	print j



gene_expression_integrated_motifs("PI_GI_Signalling_triplets.csv", "PI_GI_Signalling_quadruplets.csv", "PI_GI_Signalling_quintuplets.csv")
#gene_expression_integrated_motifs("Signalling_degree_below300_PBPF_Motifs.csv")



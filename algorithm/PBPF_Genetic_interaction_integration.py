import pandas as pd
import numpy as np
import random
import itertools
import scipy
from scipy import stats
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import pickle
from collections import defaultdict
import ast

def dictionary_index(orffile):

    orf_data=open(orffile, 'r')
    orfs=[]
    counter=0
    index={}
    for each_orf in orf_data:
        index[each_orf.strip()]=counter
        counter+=1

    inv_index = {v: k for k, v in index.iteritems()}
    print "Systematic name for YGR040W is: ", inv_index["YGR040W"] 
    return index, inv_index

dictionary_index("orf_yeast_list.csv")

def module_triplets(triplets):

	index=dictionary_index("orf_yeast_list.csv")[0]
	inv_index=dictionary_index("orf_yeast_list.csv")[1]
	print "Systematic name for YGR040W is: ", inv_index["YGR040W"]
	triplets=pd.read_csv(triplets)

	proteins_list=triplets["protein"]
	
	signalling100_triplets={}
	for protein in proteins_list:
		n3_protein=triplets.loc[triplets['protein']==protein, 'n4']
		for sides in n3_protein:
			n3_list=ast.literal_eval(sides)
			if(len(n3_list)>0):
				signalling100_triplets[protein]=(n3_list)
	
	return signalling100_triplets

module_triplets("Signalling_degree_below200_PBPF_Motifs.csv")

def genetic_interaction_of_the_triplets(genetic_interaction_matrix):

	triangle_list=module_triplets("Signalling_degree_below200_PBPF_Motifs.csv")

	
	matrix_score=np.loadtxt(genetic_interaction_matrix)
	index=dictionary_index("orf_yeast_list.csv")[0]
	m_score=np.where(np.isnan(matrix_score), 0, matrix_score)
	genetic_interaction_positive=[]
	genetic_interaction_negative=[]
	sum_matrix_score_edge=[]
	
	score_motifs=[]
	functional_triplets=[]
	for protein in triangle_list:
		triangles=triangle_list.get(protein)
		for ts in triangles:
			print ts
			side1 = ts[0]
			side2 = ts[1]
			side3 = ts[2]
			side4 = ts[3]

			#sides_list=[side1, side2, side3]
			#print sides_list
			#print side1, side2, side3
			
			matrix_score_edge1 = m_score[side1][side2]
			matrix_score_edge2 = m_score[side2][side3]
			matrix_score_edge3 = m_score[side3][side4]
			matrix_score_edge4 = m_score[side1][side4]

			#s1=matrix_score_edge1[matrix_score_edge1!=0]
			#s2=matrix_score_edge2[matrix_score_edge2!=0]
			#s3=matrix_score_edge3[matrix_score_edge3!=0]

			if(matrix_score_edge1<0.08): 
				if(matrix_score_edge2<0.08):
					if(matrix_score_edge3<0.08):

						if(matrix_score_edge4<0.08):
						#fn_triplets=[matrix_score_edge1, matrix_score_edge2, matrix_score_edge3]
							fn_triplets=[matrix_score_edge1, matrix_score_edge2, matrix_score_edge3, matrix_score_edge4]
        					#print fn_triplets
        					score_motifs.append(fn_triplets)
        					functional_triplets.append(ts)
    		print(len(functional_triplets))
    	#print "functional_triplets", functional_triplets
    	#print "score_motifs", score_motifs
	return functional_triplets

#genetic_interaction_of_the_triplets("genetic_interaction_avg_matrix2.csv")	

					#matrix_score_combined = (matrix_score_edge1, matrix_score_edge2, matrix_score_edge3)
					#print(matrix_score_combined)


import csv
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time
import seaborn as sns
import itertools
import random

orf_file="orf_yeast_list.csv"
#orf_interactions="orf_yeast_interactions.csv"
list_files=["SGA_file0.txt","SGA_file1.txt","SGA_file2.txt","SGA_file3.txt"]

num_row=6342
num_random=6

start=time.time()
def genetic_matrix(fname, orffile):

	orf_data=open(orffile, 'r')
	orfs=[]
	counter=0
	index={}
	for each_orf in orf_data:
		index[each_orf.strip()]=counter
		
		counter+=1

#	return index

#dictionary_index(orf_file)

#def genetic_matrix(fname, orffile):

	#index=dictionary_index(index)
	
	dim = len(index)
	
	matrix=np.empty([dim, dim])
	avg_matrix=np.empty([dim, dim])
	matrix[:]=np.nan
	avg_matrix[:]=np.nan
	ab_ba_different_matrix=np.empty([dim, dim])
	ab_ba_different_matrix[:]=np.nan
	repeat_mat=np.zeros((dim, dim))
	#repeat_mat[:]=np.zeros
	#print(repeat_mat)
	#return(index)
	all_source=[]
	all_target=[]
	all_source_id=[]
	all_target_id=[]
	all_score=[]
	all_pvalue=[]

	i=0
	counter=[0,0,0,0]
	for files in fname:
		raw_data_file=open(files, 'r')

		#print(files)
		counter[i]=0
		for raw_data in raw_data_file:
			if(counter[i]>0):
				each_entry=raw_data.split()
			#print(each_entry)
				presource=each_entry[0]
				pretarget=each_entry[2]
				#print(presource,'\t',pretarget)
				presource1=' _'.join(presource.split('_'))
				pretarget1=' _'.join(pretarget.split('_'))
				#print(presource1,'\t',pretarget1)
				source=presource1.split()[0]
				target=pretarget1.split()[0]
				#print(source,'\t',target)
				score=float(each_entry[5])
				p_value=float(each_entry[6])

				id_source=index[source]
				id_target=index[target]
				#print(id_source, id_target, score, p_value

#				all_source.append(source)
#				all_target.append(target)
#				all_source_id.append(id_source)
#				all_target_id.append(id_target)
#				all_score.append(score)
#				all_pvalue.append(p_value)

		

				#for ix in range(len(index)):
					#print(ix)
				iix=id_source
				iiy=id_target

				#print(iix, iiy)

				if(np.isnan(matrix[iix][iiy])):
					matrix[iix][iiy]=score
				else:
					#print(repeat_mat[iix][iiy])
					matrix[iix][iiy]+=score

				repeat_mat[iix][iiy]+=1

				#print(matrix)
	




			
			counter[i]+=1
			#print(range(counter))


		i+=1
#	all_data=pd.DataFrame({'source':all_source, 'target':all_target, 'source_id':all_source_id, 'target_id':all_target_id, 'score':all_score, 'p_value':all_pvalue})
#	all_data=all_data[['source', 'target', 'source_id', 'target_id', 'score', 'p_value']]
#	all_data.to_csv("complete_genetic_interaction_data.csv", sep='\t')

	print(matrix[0][1735])
	print(repeat_mat[0][1735])

	for aa in range(dim):
		if(repeat_mat[aa][aa]!=0):
			avg_matrix[aa][aa]=matrix[aa][aa]/repeat_mat[aa][aa]
			ab_ba_different_matrix[aa][aa]=avg_matrix[aa][aa]
		for bb in range(aa):
			if(repeat_mat[aa][bb]!=0):
			#	print(matrix[aa][bb],matrix[bb][aa])
				avg_matrix[aa][bb]=matrix[aa][bb]/repeat_mat[aa][bb]
				avg_matrix[bb][aa]=matrix[bb][aa]/repeat_mat[bb][aa]
				ab_ba_different_matrix[aa][bb]=avg_matrix[aa][bb]
				ab_ba_different_matrix[bb][aa]=avg_matrix[bb][aa]
	#max_min=np.ptp(ab_ba_different_matrix)
	#print(max_min)
	print(ab_ba_different_matrix[0][1735])
	print(avg_matrix[0][1735])
	s=[]
	t=[]
	list_ab=[]
	list_ba=[]
	

	for aa in range(dim):
		for bb in range(aa):
			if(np.isnan(matrix[aa][bb])==False and np.isnan(matrix[bb][aa])==True):
				avg_matrix[bb][aa]=avg_matrix[aa][bb]			

			elif(np.isnan(matrix[aa][bb])==True and np.isnan(matrix[bb][aa])==False):
				avg_matrix[aa][bb]=avg_matrix[bb][aa]

			elif(np.isnan(matrix[aa][bb])==False and np.isnan(matrix[bb][aa])==False):
				#print(avg_matrix[aa][bb],avg_matrix[bb][aa])
				avg_matrix[aa][bb]=(matrix[aa][bb]+matrix[bb][aa])/(repeat_mat[aa][bb]+repeat_mat[bb][aa])
				avg_matrix[bb][aa]=(matrix[aa][bb]+matrix[bb][aa])/(repeat_mat[aa][bb]+repeat_mat[bb][aa])
				
				#print(ab_ba_different_matrix[aa][bb],ab_ba_different_matrix[bb][aa])
				
				list_ab.append(ab_ba_different_matrix[aa][bb])
				list_ba.append(ab_ba_different_matrix[bb][aa])
				s.append(aa)
				t.append(bb)
				#print(matrix[aa][bb]/repeat_mat[aa][bb],matrix[bb][aa]/repeat_mat[bb][aa])
				

	print(avg_matrix[0][1735])  #ab & ba present in repeatitions
	print(avg_matrix[1735][0])
	print(avg_matrix[5470][4457])	#ab present only but in repeatition
	print(avg_matrix[4457][5470])
	print(avg_matrix[617][5137])  #ab & ba present without repeatitions
	print(avg_matrix[5137][617])
	print(avg_matrix[3350][1776])	#ab and ba present without repeatition
	print(avg_matrix[1776][3350])
	print(avg_matrix[4742][1014])	#ab present only without repeatition
	print(avg_matrix[1014][4742])

	print(avg_matrix[5517][2478])
	print(avg_matrix[2231][709])

	unique_data=pd.DataFrame({'source':s, 'target':t})
	unique_data=unique_data[['source', 'target']]
	unique_data.to_csv("Whole_genetic_interaction_with_averagescore.csv", sep='\t')

	return unique_data


	#ab_ba_list=pd.DataFrame({'source':s_id, 'target':t_id, 'ab':list_ab, 'ba':list_ba})
	#ab_ba_list=ab_ba_list[['source', 'target', 'ab', 'ba']]
	#ab_ba_list.to_csv("ab_ba_diff_score_with_genes.csv", sep='\t')


	np.savetxt("genetic_interaction_avg_matrix2.csv", avg_matrix)
	np.savetxt("genetic_interaction_ab_ba_diff_matrix3.csv", ab_ba_different_matrix)
			



	return avg_matrix, ab_ba_different_matrix, list_ab, list_ba





matrix,matrix_diff_ab_ba,ab_list,ba_list=genetic_matrix(list_files, orf_file)


print('Time = '+str((time.time()-start)/60)+'minutes')

def protein_complex(pc_infile, pc_outfile):

	index=dictionary_index(pc_infile)

	f = pd.read_csv(pc_infile)
	g = open(pc_outfile, 'w')

	n=len(f)
	i_l=[0]
	j=0.0
	for i in range(n):
		if(float(f["CID"][i])==j+1.0 and float(f["CID"][i])!=float(f["CID"][i-1])):
		#print(i)
			i_l.append(i)
			j=j+1.0	

	dict={}
	#total=0
	for i in range(len(i_l)):
		num=float(f["CID"][i_l[i]])
		#print(num)
	
		if (i==len(i_l)-1):
			end=n
		else:
			end=i_l[i+1]
	
		dict[num]=[]
		for j in range(i_l[i],end):
			dict[num].append(f["ORF"][j])
	#print len(dict[j])
	#total+=len(dict[num])

	a=[]
	b=[]


	for k in range(len(i_l)):
		#print(k)
		orfs = dict[k]
		number_orfs = len(orfs)
		if(number_orfs==4):
		#print(number_orfs)
			interactions = list(itertools.combinations(orfs, 2))
		#print interactions
			interactions_size=len(interactions)
		#print interactions_size
			 #print len(list(interactions))
			for ii in range(len(interactions)):
				c,d=interactions[ii]
				#a.append(c)
				#b.append(d)
				a.append(index[c])
				b.append(index[d])
				print(a, b)
	
	data = pd.DataFrame({'a': a, 'b': b})
	data=data[['a','b']]
    #output.write(data)
    #output.close()
	data.to_csv(g)
		#listdfs=[data_avg_file.loc[i] for i in np.split(data_avg_file.index, 40)]
			#print(listdfs)
			#for lines in listdfs:
		#ab.append(lines)
		#		scores=lines["score"]
		#		avg_scores=np.nanmean(scores)
		#		print(avg_scores)
pc_file_in = "pc_annotated.csv"
pc_file_out = "pc_annotated_outfile.csv"
#protein_complex(pc_file_in, pc_file_out) 

############################################################################################################################################
           

def orf_random(orf_outfile):

	g = open(orf_outfile, "w")

	a=[]
	b=[]
	#c=0

	for ints in range(100):
	#	a = range(num_row)
		num_row1 = 6342
		number_of_orfs = 4
		list_random=random.sample(num_row, number_of_orfs)

		interactions= list(itertools.combinations(list_random, 2))

		for genes in range(len(interactions)):
    			c,d=interactions[genes]
    			a.append(c)
    			b.append(d)
    			print(a, b)
    	
	#	c+=1
	data = pd.DataFrame({'a': a, 'b': b})
	data=data[['a','b']]
	data.to_csv(g)

out_orf = "orf_random_four_orfs.csv"
orf_random(out_orf)

def biogrid_random(bg_infile, orf_names, bg_outfile):

	bg_data=open(bg_infile, 'r')
	orf_data=open(orf_names, 'r')
	out_bg=open(bg_outfile, 'w')

	for lines in bg_data:
        	if(len(lines)>0):
                	rows=lines.split()
                	c1=rows[6]
                	c2=rows[9]
                	pre1=c1.split("gene/locuslink:")[-1]
                	pre2=c2.split("gene/locuslink:")[-1]
                	if(pre1 in list(orf_yeast_list['orfs'])):
                    		if(pre2 in list(orf_yeast_list['orfs'])):
                                #print(pre1, pre2)
                                	file.write(pre1+'\t'+pre2+'\n')

data_biogrid="BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.160.mitab.txt"
orf_yeast_list=pd.read_csv('orf_yeast_list_withheader.csv')
ppi_biogrid="biogrid_ppi.csv"
#biogrid_ppi(data_biogrid, orf_yeast_list, ppi_biogrid)


def save_random_lines_from_file(fname,fwrite,nr,nrndm):
	#fname is the file of original data file
	#fwrite is the new file where randomly selected lines from the original file fname are stored
	#nr should be the number of rows in the original file fname
	#nrndm is the number of random lines to be selected from the original file fname	
     

	x=random_array_unique_sorted(nr,nrndm) # this calls the function defined above to get an array of sorted unique random integers
	#print(x)
	f = open(fname,'r') # opens the file of original data
	g = open(fwrite,'w')	# opens file where new randomly selected data will be stored
	#SN = nx.Graph()    
	i=1 #line number in the reading file
	j=0 #index of the random array
	for line in f.readlines():
	    #print(line)
	    if(j<=(len(x)-1)):
	    	if(i==x[j]):
	    	    #print(i,j)		
	    	    g.write(line)
	            #print(line)
		    #if(x[j]==x[j+1]):
		     
        
		    

	     
	            j+=1
	             
	    i+=1
	
        f.close()
        g.close()


for files in range(100):
	
	fnam='biogrid_ppi.csv'
	frw='biogrid_pc_four'+str(files)+'.csv'
	#frw='test_sn5k'+str(files)+'.csv'
	save_random_lines_from_file(fnam,frw,num_row,num_random)
    

#pathway_list=['path1_iter.csv', 'path4_iter.csv', 'path5_iter.csv', 'path6_iter.csv', 'path8_iter.csv', 'path10_iter.csv']
#pathway_list_out=['path1_iter_score.csv', 'path4_iter_score.csv', 'path5_iter_score.csv', 'path6_iter_score.csv', 'path8_iter_score.csv', 'path10_iter_score.csv']

def average_calculation(fpcs, pcavg):


	#f=open(filepcs, 'r')
	count=0
	#for file in f:
	for file in fpcs:
		#pcs=open(file, 'r')
		pcs=pd.read_csv(file, header=None)
		
		score=[]
		a=[]
		b=[]
		#score_average=[]
		#for i in range(len(y)):
		for i in pcs:
			print(i)
			#indices=i.split(",")
			score.append(matrix[pcs[1][i]][pcs[2][i]])
			a.append(pcs[1][i])
			b.append(pcs[2][i])

			print(a, b, score)

			score_average=np.nanmean(score)
		print(score_average)

		#sns.distplot(score_average)
#plt.plot(ab_list,ba_list)
		#plt.show()
		#plt.savefig('pc4_ranodom.png')



		data = pd.DataFrame({'a': a, 'b': b, 'score': score})
		data=data[['a','b','score']]
		#data.to_csv(pathway_list_out[count])
		count+=1

		#if(score_average>0): ################ complete this to get the novel potential pcs
		#	gene_ints=

	#return score_average


#f_pcs=['pc_four'+str(files)+'.csv' for files in range(100)]
f_pcs='biogrid_pc_four1.csv'
pc_avg=['avg_score_pc_four'+str(files)+'.csv' for files in range(100)]
	#f_pcs='pc_four'+str(files)+'.csv'
	#pc_avg='avg_score_pc_four'+str(files)+'.csv'
average_calculation(f_pcs, pc_avg)





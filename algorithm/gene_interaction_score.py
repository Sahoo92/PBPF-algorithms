import pandas as pd
import csv
import numpy as np
from itertools import groupby
import itertools

genetic_interaction_data = "SGA_ExE.txt"
orf_list = "full_geni_new1.csv"
pathway_list = "list_of_genes_pathway.csv"


data_geni=open(genetic_interaction_data, 'r').readlines()
data_orfs=open(orf_list, 'r').readlines()
data_orfs1=pd.read_csv(orf_list, sep=' ')
data_pathway=pd.read_csv(pathway_list, sep=' ')
#print data_orfs1



####  INDEXING WITH DICTIONARY ####


j=0
for orfs_geni in data_orfs1:
    s = data_orfs1["Query"]
    t = data_orfs1["ID"]
    #s1=len(set(s1))
    st = s.append(t)
    st1=list(st)
    st_each_entry=(set(st))  #unique entry of each Yeast orf


    #print  [len(list(group)) for key, group in (groupby(st))]  #prints number of occurences of the interactions(A-B and B-A are different)
    ##### to be done : more than one interaction scores corresponding to the protein-interaction 
    j+=1

#print (set(s)), (set(t)), t


s1=data_orfs1[(data_orfs1['Query']=='YPR190C') & (data_orfs1['ID']=='YPL169C')]

s2=data_orfs1[(data_orfs1['Query']=='YPL169C') & (data_orfs1['ID']=='YPR190C')]


print s1
print s2


counter=0
indices=0
node_ids={}
for k in st_each_entry:
    node_ids[k]=counter
    counter+=1
    #print counter


A=np.zeros(len(st_each_entry))  
for l in st_each_entry:
    #for k_path in data_pathway["source"]:
    
    orfs = node_ids[l]
        #if(node_ids.has_key('k_path')):
    #print('the path')
    indices += 1
        #l+=1

#print orfs, st_each_entry[854]

#print node_ids #prints unique orfs Vs their assigned index


####  SCORE MATRIX  ####

n=len(node_ids)
matrix=np.zeros((n,n))
#lri=[]
repeatition={}
#keep_repeated={}
#pathway_dict={}


for i in range(len(data_orfs1)):
		ii=node_ids[data_orfs1["Query"][i]]
		jj=node_ids[data_orfs1["ID"][i]]
		repeatition[str(node_ids[data_orfs1["Query"][i]])+'-'+str(node_ids[data_orfs1["ID"][i]])]=[]

        #pathway_dict[]

for i in range(len(data_orfs1)):
		ii=node_ids[data_orfs1["Query"][i]]
		jj=node_ids[data_orfs1["ID"][i]]


		if(matrix[ii][jj]!=0):
			#print(ii,jj)
			#print(repeatition[str(ii)+'-'+str(jj)])		
			
			repeatition[str(ii)+'-'+str(jj)].append(data_orfs1["name"][i])

			if(len(repeatition[str(ii)+'-'+str(jj)])==1):
				repeatition[str(ii)+'-'+str(jj)].append(matrix[ii][jj])
			
			#print(repeatition[str(ii)+'-'+str(jj)])
			#repeatition[data_orfs1["Query"][i]]

		matrix[ii][jj]= data_orfs1["name"][i]

		

                #repeatition[data_orfs1["Query"][i]+'-'+data_orfs1["ID"][i]]=[]
                #repeatition[data_orfs1["Query"][i]+'-'+data_orfs1["ID"][i]].append(data_orfs1["name"][i])



		#	matrix[ii][jj]= data_orfs1["name"][i]
        


#print repeatition

#print(data_orfs1)

print matrix[237][353]
np.savetxt("geni_matrix.csv", matrix, delimiter=" ")
#print(matrix[5][600], matrix[600][5])
#print(lri)
#print(repeatition)
nz=np.ndarray.nonzero(matrix)
#high_score=np.ndarray.argmax([matrix, 1]

#lll=0

##### A-B Vs B-A SCORE ######

d1=open("ab_geni_int.csv", "w")
d2=open("ba_geni_int.csv", "w")

for i in range(len(data_orfs1)):
        ii=node_ids[data_orfs1["Query"][i]]
        jj=node_ids[data_orfs1["ID"][i]]


        if(repeatition.has_key(str(ii)+'-'+str(jj))):
            if(repeatition.has_key(str(jj)+'-'+str(ii))):
                if(len(repeatition[str(ii)+'-'+str(jj)])>0):
                    if(len(repeatition[str(jj)+'-'+str(ii)])>0):

            
                        long_list1=repeatition[str(ii)+'-'+str(jj)]
                        long_list2=repeatition[str(jj)+'-'+str(ii)]

                        summation1=sum(long_list1)
                        summation2=sum(long_list2)

                        average1=np.mean(long_list1)
                        average2=np.mean(long_list2)


                        #print average1, average2

                        d1.write(str(average1)+"\n")
                        d2.write(str(average2)+"\n")



                        #with open('ab_ba_geni_int.csv', 'w') as f:
                        #    writer = csv.writer(f, delimiter=' ')
                        #    writer.writerows(zip(float(average1), float(average2)))

            #lll+=1


###### Set of genes ----- Pathway ######

d11=open("ab_pathway.csv", "w")
d22=open("ba_pathway.csv", "w")

ll=0
for genes in data_pathway:

    #print data_pathway["source"]

    pathway_ints = list(itertools.combinations(data_pathway["source"], 2))

    #ints_all = pathway_ints[1].split()[0]
    for i in range(len(pathway_ints)):
        #pathway_ints[i]
        a,b=pathway_ints[i]
        m=node_ids[a]
        n=node_ids[b]
        if(repeatition.has_key(str(m)+"-"+str(n))):
            if(repeatition.has_key(str(n)+"-"+str(m))):
                if(len(repeatition[str(m)+'-'+str(n)])>0):
                    if(len(repeatition[str(n)+'-'+str(m)])>0):

                        list1=repeatition[str(m)+"-"+str(n)]
                        list2=repeatition[str(n)+"-"+str(m)]

                        summation11=sum(list1)
                        summation22=sum(list2)

                        #a1=summation11/float(len(list1))
                        #a2=summation22/float(len(list2))

                        average11=np.mean(list1)
                        average22=np.mean(list2)
      
            #list1=[int(sum(repeatition[str(m)+"-"+str(n)])) for sum(repeatition[str(m)+"-"+str(n)]) in list1]
                
                        #print a1, average11
                        #print a2, average22
                        #print average11, average22

                        d11.write(str(average11)+"\n")
                        d22.write(str(average22)+"\n")

                        ll+=1

                
#####  RAW DATA CLEANING  #####  

#data_protein-complex_orfs=open()

data1=pd.read_csv("ab_geni_int.csv", sep=' ')
data2=pd.read_csv("ba_geni_int.csv", sep=' ')

merged1=pd.merge(data1, data2, left_index=True, right_index=True)
data12=merged1.to_csv("abVsba_geni_int.csv", index=False)

data3=pd.read_csv("ab_pathway.csv", sep=' ')
#print data11
data4=pd.read_csv("ba_pathway.csv", sep=' ')

merged2=pd.merge(data3, data4, left_index=True, right_index=True)
#print merged2
data34=merged2.to_csv("abVsba_pathway.csv", index=False)


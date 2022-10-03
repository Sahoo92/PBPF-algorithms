import pandas as pd  
import numpy as np
import itertools

data = pd.read_csv("pc_annotated.csv")

#print(data[["CID"]])
#print(float(data["CID"][22]))

n=len(data)

#print(n)

i_l=[0]
#i_l.append(1)
#print(i_l)


j=0.0
for i in range(n):
	if(float(data["CID"][i])==j+1.0 and float(data["CID"][i])!=float(data["CID"][i-1])):
		#print(i)
		i_l.append(i)
		j=j+1.0	

#print(len(i_l))		

dict={}
#total=0
for i in range(len(i_l)):
	num=float(data["CID"][i_l[i]])
	#print(num)
	
	if (i==len(i_l)-1):
		end=n
	else:
		end=i_l[i+1]
	
	dict[num]=[]
	for j in range(i_l[i],end):
		dict[num].append(data["ORF"][j])
        #print len(dict[j])
	#total+=len(dict[num])

#print(total)
f1=open("size_of_protein_complexes.csv", "w")
f2=open("orfs_of_protein_complexes.csv", "w")
f3=open("number_of_interactions_in_protein_complexes.csv", "w")
for k in range(len(i_l)):

    
    orfs = dict[k]
    number_orfs = len(orfs)
    interactions = itertools.combinations(orfs, 2)
    interactions_size=len(list(interactions))
    #print len(list(interactions))

    f1.write(str(number_orfs)+"\n")
    f2.write(str(orfs)+"\n")
    f3.write(str(interactions_size)+"\n")





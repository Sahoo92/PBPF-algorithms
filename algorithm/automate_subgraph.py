import csv
import networkx as nx
import pandas as pd
import numpy as np
import random
import sys
#import time


yeast_regnet = 'yeast_reg1b.csv'


def load_network(inputfile):
    """
     function to load inputfile into graph
     """


    G=nx.Graph()

    with open('yeast_reg1b.csv', 'rb') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            source = row[0]
            target = row[1]
            words1 = source.split()
            words2 = target.split()
            words1.sort()
            words2.sort()

            for word1 in words1:
                for word2 in words2:

     
                  #print word1, word2



                 G.add_edge(word1, word2)

                 #deg1 = G.degree(ABF1)
                 #reg1 = nx.degree(G, ABF1)

                 #print len(reg1) 

    return(G) 
 
G_RN = load_network(yeast_regnet)

names = G_RN.nodes()
reg1 = nx.degree(G_RN)

#print reg1


counter = 0
indices=0
node_ids = {}

for i in names:
    node_ids[i] = counter
    counter += 1
    
print len(names), counter

H = nx.Graph()


A=np.zeros(len(G_RN.edges()))   
B=np.zeros(len(G_RN.edges()))   
for i in G_RN.edges():
    
    source_a = i[0]
    target_b = i[1]

    source_a_indices = node_ids[source_a]
    target_b_indices = node_ids[target_b]

    

    H.add_edge(source_a_indices, target_b_indices)



   
    A[indices]=source_a_indices       
    B[indices]=target_b_indices   
    #empty = np.zeros(shape=(6726,2))
    #for i in range(6726):
    	     #empty[i,0]=a_new[i]
    	     #empty[i,1]=b_new[i]

    	     #print empty[i,0]
    #print A[indices] 
    #for indices in A:
#sg_data=random.sample(A[indices], 30)
         #print indices

         #S_G = nx.subgraph(H,sg_data)         
    #sg_data=random.sample(A, 30)
    
    indices += 1
    
    with open("stringsort_sid.py") as fp:
       for j, line in enumerate(fp):
        if "\xe2" in line:
            print j, repr(line)
reg2 = nx.degree(H)
#print reg2

print "Degree of LYS1 in first network", reg1['LYS1']
ii = node_ids['LYS1']
print "Degree of LYS1 in second network", reg2[ii]


     


f = open("yeastreg_data_outputa.csv",'w')

for i in xrange(len(A)):
    f.write("{} {}\n".format(A[i].astype(int), B[i].astype(int)))   

f.close()

def load_subnetwork(n):
     
     S_G = nx.Graph()
 

#loop = 0
     sg_source=random.sample(A, 30)
     sg_target=random.sample(B, 30)

     #S_G.add_edge(sg_source, sg_target)
     

#loop+=1
     
     
     return(S_G)

    

#names_sg = S_G.nodes()
#print names_sg
#print sg_data


    #print(time.time(), time.clock())
#def load_subnetwork(inputfile):    

     #S_G = nx.Graph()
         
     #with open('yeastreg_data_outputa.csv', 'rb') as sg:
         #sub_graph = csv.reader(sg, delimiter=';')
         #for row in sub_graph:
             #source_subgraph = row[0]
             #target_subgraph = row[1]
            
             #G.add_edge(source_subgraph, target_subgraph)

     #return(S_G)

"""Collect command-line options in a dictionary"""

def optparse(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

if __name__ == '__main__':
    from sys import argv
    myargs = getopts(argv)
    if '-i' in myargs:  # Example usage.
        print(myargs['-i ./fanmod_command_line_linux 3 100000 1 yeastreg_data_outputa.csv 1 0 0 2 0 0 0 1000 3 3 yreg_ms3 0 0'])
    print(myargs)  
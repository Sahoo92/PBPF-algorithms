import numpy as np
import networkx as nx
import subprocess
from subprocess import  Popen, PIPE, STDOUT
import matplotlib.pyplot as plt

num_row=21157
num_random=num_row


#To be able to use the below functions in a file, keep this file in the same folder and write the following lines in the top of the file
#import test #(or whatever is the name of this file)
#import test.save_random_lines_from_file

#then call the function as 
#save_random_lines_from_file(name of original file, name of new file, num of lines in the original file, no. of random lines to be taken)

def random_array_unique_sorted(nr,nrndm):



	xr=np.zeros(nrndm)
	for i in range(nrndm):
		if(i==0):
			xr[i]=np.random.random_integers(1,nr,1)
		if(i!=0):
			z=1
			while(z!=0):
				xr[i]=np.random.random_integers(1,nr,1)
				z=0
				for j in range(i):
					if(xr[j]==xr[i]):
						z+=1			
	xr=np.sort(xr)
	return(xr)

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


for files in range(20):
	#fnam='yeastreg_data_outputa.csv'
	fnam='path4_index1.csv'
	frw='path4'+str(files)+'.csv'
	#frw='test_sn5k'+str(files)+'.csv'
	save_random_lines_from_file(fnam,frw,num_row,num_random)
    


	#cmd = './fanmod_command_line_linux 5 100000 1 test_sn5k'+str(files)+'.csv 1 0 0 2 0 0 0 1000 3 3 yreg_ms55k'+str(files)+'.csv 0 0'
        #p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

	cmd = './fanmod_command_line_linux 3 100000 1 path4'+str(files)+'.csv 1 0 0 2 0 0 0 1000 3 3 yeast_path4'+str(files)+'.csv 0 0'
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    
	output = p.stdout.read()
	print output
##_command_line_linux

def read_fanmod_files(fan_file,motif_file,row,column):
    
    
	f=open(fan_file, 'r')
	g=open(motif_file, 'w')
	i=1
	#j=1
	for line in f:
		parts = line.split() # split line into parts
		if(i==row):
		    #print parts
		    #if(j==parts[column]):

		    	#g.write(line)
		    	#print j
			print(parts[column-1])
			
			g.write(parts[column-1])
			#print g
        #print(i)
		i+=1
    	f.close()
        g.close()



for data_files in range(20):
	fan_file='yeast_path4'+str(data_files)+'.csv'
	#print fan_file
	read_fanmod_files(fan_file,'data_motif.csv',18,1)

    
    #read_fanmod_files('yreg_ms3'+str(files)+'.csv','data_motif.csv',18,1)	
#plt.plot(parts[column-1])
#plt.axis([0,10,0,60])
#plt.show()



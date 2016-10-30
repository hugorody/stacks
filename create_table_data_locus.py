#!/usr/bin/python
#usage: python create_table_large_locus.py file.tags.tsv batch_1.catalog.tags.tsv
#CREATE TABLE WITH DATA FOR EACH LOCUS
#created by H.Rody Oct 30, 2016
'''
Columns of the table
idlocus = individual_locus
sequenc = locus consensus sequence
logqual = Log likelyhood of the locus
blackli = whether locus is blacklisted (1) or not (0)
ncopies = number of stacks for the locus
nsegind = number of individuals in the population sharing the locus
'''

import sys
import time


filetags = sys.argv[1]
filebatch = sys.argv[2]

try:
    set_tags = open(filetags)
except IOError:
    print ("File TAGS doesn't exist!")


for i in range(1):
    line=set_tags.next().split("\t")
    individual = line[1]

print "Creating catalog1 and catalog2...",time.ctime()

catalog = {}
catalog2 = {}       #number of stacks

for line in set_tags:
	line = line.split("\t")
	idlocus = line[1]+"_"+line[2]
	
	if "consensus" in line:
		
		sequenc = line[9]
		logqual = line[13].rstrip()
		blackli = line[11]
				
		mylist = []
		mylist.append(sequenc)
		mylist.append(logqual)
		mylist.append(blackli)
		
		catalog[idlocus] = mylist

	if "primary" in line:
		classinum = line[6]+line[7]

		if idlocus in catalog2:											#catalog2
			catalog2[idlocus].append(classinum)
		else:
			catalog2[idlocus] = [classinum]

set_tags.close()

########################################################################
set_batch = set(line.strip() for line in open (filebatch, 'r'))

print "Creating catalog3...",time.ctime()

catalog3 = {}      #<id> <individuals>

for i in set_batch:
	i = i.split("\t")

	seg = i[8].split(",")

	lista = []
	
	for j in seg:
		j = j.split("_")
		k = j[0]

		lista += k.splitlines()

	for z in seg:
		catalog3[z] = len(set(lista))  #add id plus number of individuals

########################################################################

print "Merging catalogs and creating final table...",time.ctime()

output1 = open(individual+".txt", "w")
output1.write("#ID_LOCUS SEQUENCE LOG_LIKELIHOOD NUMBER_STACKS NUMBER_IND_SHARING\n");   #write header in output1

for x in catalog.keys():
	if x in catalog2.keys():
		if x in catalog3.keys():
			#print x," ".join(catalog[x]),len(set(catalog2[x])),catalog3[x]
			output1.write(str(x)+" "+" ".join(catalog[x])+" "+str(len(set(catalog2[x])))+str(catalog3[x])+"\n");  #write data output1

output1.close()
set_batch.close()

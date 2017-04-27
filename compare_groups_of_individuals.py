#!/usr/bin/python
#usage: 

import sys
import itertools

file1 = sys.argv[1] #a list containing sampleIDs separated by line break
file2 = sys.argv[2] #batch_1.catalog.tags.tsv


set_1 = list(line.strip() for line in open (file1, 'r'))
set_2 = list(line.strip() for line in open (file2, 'r'))

nameoutput = file1.replace(".lst", "")

output0 = open(nameoutput+"_table.txt", "w")  #OUTPUT 0: TABLE
output0.write("Desc Nshared_loci\n");

#generate group catalog
output1 = open(nameoutput+"_group.txt", "w")  #OUTPUT 1 GROUP CATALOG
output1.write("CatalogID SampleID_LocusID\n");

sizegroup = []
for i in set_2:
	if "#" not in i:
		i = i.split("\t")
		
		catalogid = i[2] #id number of matches in catalog
		listsamplelocus = i[8].split(",") #sampleID_locusID list of the matches
		
		dictsampleid = {} #dictionary containing the sampleID as keys, and locusID as values
		
		for j in listsamplelocus:
			j = j.split("_") 
			
			sampleid1 = j[0]
			locusid1 = j[1]
			
			dictsampleid[sampleid1] = locusid1 #add sampleidual id to listindiv
		
		newdict = {}  #dictionary containing keys and values for correspondent samplesIDs from mylist file
		
		for i in set_1:
			if i in dictsampleid:
				newdict[i] = dictsampleid[i]

		
		if int(len(newdict)) == int(len(set_1)):
			sizegroup.append(catalogid)
			
#			print catalogid+" "+",".join(["_".join([key, str(val)]) for key, val in newdict.items()])
			output1.write(catalogid+" "+",".join(["_".join([key, str(val)]) for key, val in newdict.items()])+"\n");


output0.write("group "+str(len(sizegroup))+"\n");

#generate catalog for pair combinations among group members

pairs = [] #all pair combination of individuals from file1
for pair in itertools.combinations(set_1,2):
	pairs.append(pair)

for sampleid in pairs:
	id1 = sampleid[0]
	id2 = sampleid[1]
	
	sizepair = []
	
	output2 = open(nameoutput+"_"+id1+"_"+id2+".txt", "w")  #OUTPUT 2 INDIVIDUAL PAIR CATALOG
	output2.write("CatalogID SampleID_LocusID\n");
	
	for i in set_2:
		if "#" not in i:
			i = i.split("\t")
			
			catalogid = i[2] #id number of matches in catalog
			listsamplelocus = i[8].split(",") #sampleID_locusID list of the matches
			
			dictsampleid = {} #dictionary containing the sampleID as keys, and locusID as values
			
			for j in listsamplelocus:
				j = j.split("_") 
				
				sampleid1 = j[0]
				locusid1 = j[1]
				
				dictsampleid[sampleid1] = locusid1 #add sampleidual id to listindiv
				
			
			if id1 in dictsampleid and id2 in dictsampleid:
				
#				print catalogid,id1+"_"+dictsampleid[id1],id2+"_"+dictsampleid[id2]
				output2.write(catalogid+" "+id1+"_"+dictsampleid[id1]+","+id2+"_"+dictsampleid[id2]+"\n");
				sizepair.append(catalogid)
	
	output0.write(id1+","+id2+" "+str(len(sizepair))+"\n");

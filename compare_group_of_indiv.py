#!/usr/bin/python
#usage: python compare_group_of_individuals.py list.lst batch_1.catalog.tags.tsv int
#sumarize the interactions of a subpopulation of individuals of a larger population/group of individuals to help understand if stacks analysis are making sense
#creates four outputs: output0 _table.txt has two columns Desc (Group, for all individuals; and Indiv1,Indiv2 for pairs of individuals) and Nshared_loci (). output1 _group.txt has two columns CatalogID (Catalog ID of the match) and SampleID_LocusID (Individual ID plus its Locus ID combined with _ and separated by comma). 
#output2 has the same columns of output1 but now showing the locus information that is shared by all possible individual pair combinations
#outputfasta is fasta sequences for loci shared among groups

import sys
import itertools

file1 = sys.argv[1] #a list containing sampleIDs separated by line break
file2 = sys.argv[2] #batch_1.catalog.tags.tsv
nindi = sys.argv[3] #minimum number (int) of individuals to consider when calculating the catalog

set_1 = list(line.strip() for line in open (file1, 'r'))
set_2 = list(line.strip() for line in open (file2, 'r'))

outputdir = sys.argv[4]

output0 = open(outputdir+"compare_table.grp", "w")  #OUTPUT 0: TABLE
output0.write("Desc Nshared_loci\n");

#generate group catalog
output1 = open(outputdir+"compare_group.grp", "w")  #OUTPUT 1 GROUP CATALOG
output1.write("CatalogID SampleID_LocusID\n");

#generate fasta for group catalog
outputfasta = open(outputdir+"sequences.fasta","w")


wholecatalog = {}

sizegroup = []
for i in set_2:
	if "#" not in i:
		i = i.split("\t")
		
		catalogid = i[2] #id number of matches in catalog
		sequence = i[9]
		listsamplelocus = i[8].split(",") #sampleID_locusID list of the matches
		
		dictsampleid = {} #dictionary containing the sampleID as keys, and locusID as values
		
		wholecatalog[catalogid] = sequence
				
		for j in listsamplelocus:
			j = j.split("_") 
			
			sampleid1 = j[0]
			locusid1 = j[1]
			
			dictsampleid[sampleid1] = locusid1 #add sampleidual id to listindiv
		
		newdict = {}  #dictionary containing keys and values for correspondent samplesIDs from mylist file
		
		for i in set_1:
			if i in dictsampleid:
				newdict[i] = dictsampleid[i]

		
#		if int(len(newdict)) == int(len(set_1)):  #filter based in the percentage of individuals
		if int(len(newdict)) >= int(nindi):  #filter based in the number of individuals nindi variable
			sizegroup.append(catalogid)
			
#			print catalogid+" "+",".join(["_".join([key, str(val)]) for key, val in newdict.items()])
			output1.write(catalogid+" "+",".join(["_".join([key, str(val)]) for key, val in newdict.items()])+"\n");
			outputfasta.write(">"+catalogid+"\n"+str(wholecatalog[catalogid])+"\n");
			

print "group "+str(len(sizegroup))
output0.write("group "+str(len(sizegroup))+"\n");

#generate catalog for pair combinations among group members

pairs = [] #all pair combination of individuals from file1
for pair in itertools.combinations(set_1,2):
	pairs.append(pair)

for sampleid in pairs:
	id1 = sampleid[0]
	id2 = sampleid[1]
	
	sizepair = []
	
	output2 = open(outputdir+id1+"_"+id2+".grp", "w")  #OUTPUT 2 INDIVIDUAL PAIR CATALOG
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
	
	print id1+","+id2+" "+str(len(sizepair))
	output0.write(id1+","+id2+" "+str(len(sizepair))+"\n");
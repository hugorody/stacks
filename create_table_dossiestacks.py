#!/usr/bin/python
#read stacks outputs and calculate some indices for comparison among different runs
#usage: python create_table_dossiestacks.py file.tags.tsv file.alleles.tsv batch_1.catalog.tags.tsv

import sys

filetags = sys.argv[1]
filealle = sys.argv[2]
filebatc = sys.argv[3]

try:
    set_tags = open(filetags)
except IOError:
    print ("File TAGS doesn't exist!")

set_alle = set(line.strip() for line in open (filealle, 'r'))
set_batc = set(line.strip() for line in open (filebatc, 'r'))



###################################################### FILE 1

#Total number of locus or consensus stacks
number_consensus = 0

#Number of blacklisted tags
blacklisted = []

#Used for Number_of_loci_with_copies_>_14
catalog2 = dict() 

for line in set_tags:
	line = line.split("\t")
	log_likelihood = line[13]

	if "consensus" in line:
		if float(log_likelihood) >= -50.0:    #log likelihood cutoff for loci greater than -50
			number_consensus += 1
		
		if int(line[11]) == 1:
			blacklisted += str(line[1]+"_"+line[2]).splitlines()


	
	if "primary" in line:# or "secondary" in line:
		idstack = line[1]+"_"+line[2]
		classinum = line[6]+line[7]
		
		if idstack in catalog2:
			catalog2[idstack].append(classinum)
		else:
			catalog2[idstack] = [classinum]

#Number_of_loci_with_copies_>_14
number_of_copies = 0

for i in catalog2.items():
	idstack = i[0]
	ncopies = len(set(i[1]))
	
	if int(ncopies) > 14:
		number_of_copies += 1
	
	
	
	
print "Total_number_of_predicted_loci",number_consensus
print "Number_of_blacklisted_loci",len(set(blacklisted))
print "Number_of_loci_with_copies_>_14",number_of_copies

###################################################### FILE 2

#Number of locus with predicted alleles in the lineage
locus_alle = []

for line1 in set_alle:
	line1 = line1.split("\t")
	
	locus_alle += line1[2].splitlines()

	

print "Number_of_loci_with_predicted_alleles",len(set(locus_alle))
print "Mean_number_of_predicted_alleles_per_locus",float(float(len(set_alle))/float(len(set(locus_alle))))


###################################################### FILE 3

#Number of segretating locus
segregating_locus_high = 0
segregating_locus_low = 0
privative_loci = 0


for line2 in set_batc:
	line2 = line2.split("\t")
	
	distribution = len(line2[8].split(","))
	
	if distribution >= 166:
		segregating_locus_high += 1
		
	if distribution >= 14:
		segregating_locus_low += 1
	
	if distribution == 1:
		privative_loci += 1
		
		
print "Number_of_segregating_loci_(>=_166)",segregating_locus_high
print "Number_of_segregating_loci_(>=_10)",segregating_locus_low
print "Number_of_privative_loci",privative_loci

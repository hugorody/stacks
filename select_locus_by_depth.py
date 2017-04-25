#!/usr/bin/python
#usage: python select_locus_by_depth.py file.matches.tsv file.models.tsv loglikelihood_value depth_value

import sys
import numpy

file1 = sys.argv[1] #.matches.tsv: Matches to the catalog
file2 = sys.argv[2] #models.tsv
file3 = sys.argv[3] #alleles.tsv
file4 = sys.argv[4] #batch_1.catalog.tags.tsv
var1 = sys.argv[5]  #minimum log likelihood cutoff
var2 = sys.argv[6]  #minimum stack depth

logllcut = float(var1)   # if a locus is a perfect homozygote, with no sequencing errors, each column will produce a likelihood of 1. That is, a 100% likelihood of a homozygote. Ln(1) = 0, so the best score a locus can have in log likelihood world is 0, and as the locus gets worse the likelihood gets more negative.
depthlcut = int(var2)

try:
	set_1 = open(file1)
	set_2 = open(file2)
	set_3 = open(file3)
	set_4 = open(file4)
except IOError:
	print ("Input files doesn't exist!")


selectedconsensus = {}  #dictionary containing SampleID_LocusID as keys, and Stack depth as values

for i in set_1:   #.matches.tsv: Matches to the catalog
	i = i.split("\t")
	
	if "consensus" in i:
		depthl = int(i[6]) #stack depth. number or reads contained in the locus that matched to the catalog.
		logll = float(i[7])  #Log likelihood of the matching locus.
		
		if depthl >= depthlcut and logll >= logllcut: #locus filtering. Pass only locus with depth >= cutoff (var2), and log likelihood greater than cutoff (var1)
			
			selectedconsensus[i[3]+"_"+i[4]] = depthl


alleles = []

for j in set_3:
	if "#" not in j:
		j = j.split("\t")
	
		samplelocusid = j[1]+"_"+j[2]
		
		if samplelocusid not in alleles:
			alleles.append(samplelocusid)



nlocusamb = int(len(set(alleles)))  #number of stacks/locus with SNPs
nlocus = int(len(selectedconsensus)) #total number of stacks/locus 
percent = (float(nlocusamb) * 100.0) / nlocus #percentage of locus presenting SNPs

catalog = {}

for u in set_4: #batch_1.catalog.tags.tsv
	if "#" not in u:
		u = u.split("\t")
		key = u[2]
		line = u[8]
		numtimes = int(len(line.split(","))) #number of individuals the locus is shared
		linelist = line.split(",") #list of sampleID_locusID being shared
		
		for g in linelist:
			catalog[g] = numtimes
	
sharetimes = [] #list containing the number of times each locus presenting SNPs is shared across the individuals

for l in alleles: #for each of locus presenting SNPS
	if l in catalog: #if the locus is in the catalog
		sharetimes.append(catalog[l]) #add to a list the number of individuals share this locus
	
average = numpy.mean(sharetimes) #for the pool of locus presenting SNPs, the global average number of individuals sharing those locus

print nlocus,nlocusamb,percent,average

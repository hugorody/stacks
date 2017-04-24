#!/usr/bin/python
#usage: python select_locus_by_depth.py file.matches.tsv file.models.tsv loglikelihood_value depth_value

import sys

file1 = sys.argv[1] #matches.tsv: Matches to the catalog
file2 = sys.argv[2] #models.tsv
var1 = sys.argv[3]  #minimum log likelihood cutoff
var2 = sys.argv[4]  #minimum stack depth

logllcut = float(var1)   # if a locus is a perfect homozygote, with no sequencing errors, each column will produce a likelihood of 1. That is, a 100% likelihood of a homozygote. Ln(1) = 0, so the best score a locus can have in log likelihood world is 0, and as the locus gets worse the likelihood gets more negative.
depthlcut = int(var2)

try:
	set_1 = open(file1)
	set_2 = open(file2)
except IOError:
    print ("Input files doesn't exist!")


selectedconsensus = {}  #dictionary containing SampleID_LocusID as keys, and Stack depth as values

for i in set_1:   #matches.tsv: Matches to the catalog
	i = i.split("\t")
	
	if "consensus" in i:
		depthl = int(i[6]) #stack depth. number or reads contained in the locus that matched to the catalog.
		logll = float(i[7])  #Log likelihood of the matching locus.
		
		if depthl >= depthlcut and logll >= logllcut: #locus filtering. Pass only locus with depth >= cutoff (var2), and log likelihood greater than cutoff (var1)
			
			selectedconsensus[i[3]+"_"+i[4]] = depthl


print len(selectedconsensus)

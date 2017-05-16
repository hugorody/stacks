#!/usr/bin/python
#verify the frequency of predicted SNPs in each base of all loci of individuals and create a scatter plot
#input: *.tags.tsv (only for individuals files, not batch)
#usage: python script.py *.tags.tsv

import sys
import matplotlib.pyplot as plt
import numpy as np

file1 = sys.argv[1]  #tags.tsv
figname = file1.replace(".tags.tsv",".png")

set1 = list(line.strip() for line in open (file1, 'r'))

models = dict()  #contains a concise record of all the model calls for each nucleotide in the locus ("O" for homozygous, "E" for heterozygous, and "U" for unknown)
for i in set1:
	if "model" in i:
		i = i.split("\t")
		idsample = i[1]+"_"+i[2]
		snpseq = i[9]
		
		models[idsample] = snpseq


chart = dict()  #each nucl position of all individual loci as key and the SNP frequency as values 
for i in models.items():
	if "E" in i[1]:
		positions = [pos for pos, char in enumerate(i[1]) if char == "E"]  #list containing positions where the index E occurs
		for ii in positions:
			if ii in chart:
				soma = int(chart[ii]) + 1
				chart[ii] = soma
			else:
				chart[ii] = 1
				
mybases = []
myhits = []
		
for i in chart.items():
	mybases.insert(0, i[0]) #create ordered list for nucl positions. will be used as x axis
	myhits.insert(0, i[1]) #create ordered list for SNP ferquency. will be used as y axis

numloci = len(models) #total number of loci

plt.plot(mybases, myhits)
plt.xlabel('Bases')
plt.ylabel('Frequency SNPs')
plt.savefig(figname)
#plt.show()
print "end",figname,numloci

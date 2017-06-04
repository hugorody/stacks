#!/usr/bin/python
#usage: python select_locus_by_depth.py file.matches.tsv file.models.tsv loglikelihood_value depth_value
#Calculates parameters to summarize Stacks results for a specific individual, based on a filter of minimum stack depth and a loglikelihood cutoff (normally used as greater than -10)
#Detailed description is given at the end of script, in the Variables calc section 

import sys
import numpy

file1 = sys.argv[1] #matches.tsv: Matches to the catalog
file2 = sys.argv[2] #tags.tsv
file3 = sys.argv[3] #alleles.tsv
file4 = sys.argv[4] #batch_1.catalog.tags.tsv
var1 = sys.argv[5]  #minimum log likelihood cutoff
var2 = sys.argv[6]  #minimum stack depth

logllcut = float(var1)   # if a locus is a perfect homozygote, with no sequencing errors, each column will produce a likelihood of 1. That is, a 100% likelihood of a homozygote. Ln(1) = 0, so the best score a locus can have in log likelihood world is 0, and as the locus gets worse the likelihood gets more negative. Normally is used greater than -10.
depthlcut = int(var2)

try:
	set_1 = open(file1)
	set_3 = open(file3)
	set_4 = open(file4)
except IOError:
	print ("Input files doesn't exist!")

#File 1: matches.tsv
locusmatches = []  #list containing SampleID_LocusID of locus that is present on the matches list (matched to catalog) and passed cutoff filter

for i in set_1:   #matches.tsv: Matches to the catalog
	i = i.split("\t")
	
	if "consensus" in i:
		locusmatch = i[3]+"_"+i[4]
		depthlocus = int(i[6])
		loglocus = float(i[7])
		
		if depthlocus >= depthlcut and loglocus >= logllcut:
			locusmatches.append(locusmatch)

set_1.close()

#File 3: alleles.tsv
alleles = [] #list containing all locus containing predicted alleles (SNPs)

for j in set_3:  #alleles.tsv
	if "#" not in j:
		j = j.split("\t")
	
		alleleid = j[1]+"_"+j[2] #merge sampleID_locusID
		
		if alleleid not in alleles:
			alleles.append(alleleid)

set_3.close()

#File 4: batch_1.catalog.tags.tsv
catalog = {} #dictionary containing sampleID_locusID as keys, and the number of times the locus is shared among individuals as values

for u in set_4: #batch_1.catalog.tags.tsv
	if "#" not in u:
		u = u.split("\t")
		key = u[2]
		line = u[8]
		numtimes = int(len(line.split(","))) #number of individuals the locus is shared
		linelist = line.split(",") #list of sampleID_locusID being shared
		
		for g in linelist:
			catalog[g] = numtimes
					
set_4.close()

#File 2: tags.tsv

set_2 = set(line.strip() for line in open (file2, 'r'))

consensusqual = {}  #sampleID_locusID loglikelihood
consensusdepth = {} #sampleID_locusID depth_reads
consensusselec = [] #sampleID_locusID list where loglikelihood and depth passed the cutoff

for k in set_2:
	if "consensus" in k:
		k = k.split("\t")
		samlocus = k[1]+"_"+k[2]
		locusqual = k[13]
		consensusqual[samlocus] = locusqual
		consensusdepth[samlocus] = 0

for k in set_2:
	if "primary" in k:
		k = k.split("\t")
		samlocus = k[1]+"_"+k[2]
		
		if samlocus in consensusdepth:
			consensusdepth[samlocus] += 1

for k in consensusqual.items():
	locu = k[0]
	qual = float(k[1])		
	dept = int(consensusdepth[locu])
	
	if qual >= logllcut and dept >= depthlcut:
		consensusselec.append(locu)

del set_2

#### Lists

#List
globaldepth = [] #containing the depth number for each predicted locus that passed cutoff filter

for i in consensusdepth.items():
	if i [0] in consensusselec:
		globaldepth.append(i[1])

#List 4
snpdepth = [] #containing the depth number for each predicted locus presenting SNPs 

for i in alleles:
	if i in consensusselec:
		snpdepth.append(consensusdepth[i])

#List 2
sharetimes = [] #list containing the number of times each locus presenting SNPs is shared across the individuals

for l in alleles: #for each of locus presenting SNPS
	if l in catalog: #if the locus is in the catalog
		sharetimes.append(catalog[l]) #add to a list the number of individuals share this locus

#Variables calc
totalpredictedlocus = len(consensusqual) #total number of predicted locus
totalselectedlocus = len(consensusselec) #total number of predicted locus that passed cutoff filter
nlocusmatches = len(locusmatches) #total number of stacks/locus that matched to catalog
nlocusamb = len(set(alleles))  #number of stacks/locus with SNPs
avedepthamblocus = numpy.mean(snpdepth) #the number of locus presenting SNPs with depth greater than cutoff 
percentamblocus = (float(nlocusamb) * 100.0) / totalselectedlocus #percentage of locus presenting SNPs
avedepthglobal = numpy.mean(globaldepth) #the global average depth number for the pool of locus with depth greater than cutoff
avesharesnp = numpy.mean(sharetimes) #for the pool of locus presenting SNPs, the global average number of individuals sharing those locus

print totalpredictedlocus,totalselectedlocus,avedepthglobal,nlocusmatches,nlocusamb,avedepthamblocus,percentamblocus,avesharesnp
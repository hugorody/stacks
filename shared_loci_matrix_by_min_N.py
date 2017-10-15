#!/usr/bin/python
#This script is the part 1st of "From Stacks to PCA" pipeline
#Uses batch_1.catalog.tags.tsv output file from Stacks to generate a csv space delimited matrix file where the lines are the individuals and columns the loci.
#If an individual shares the respective loci, the value is 1 whereas for the opposite the value is 0.

#As a second input this script also takes a tab delimited population file <isolatecode><tab><lineage><tab><isolate><tab><idstacks>
#Finally, a minimum number of samples in a locus must be set as parameter

#Input: 2 files and 1 number int
#------ 1. Stacks file "batch_1.catalog.tags.tsv"
#------ 2. Population file with 4 columns tab delimited
#------------- 1 column: Individual Code, for example MY_CODE003
#------------- 2 column: Individual lineage, if any it can be set as 1 for all individuals
#------------- 3 column: Individual file name, for example 003
#------------- 4 column: Individual Stacks ID
#------ 3. Int value for a minimum number of individuals sharing the locus
#Outputs:
#------ 1. CSV file named samplelocimatrix.csv

import sys
file1 = sys.argv[1]
file2 = sys.argv[2]
nsamples = sys.argv[3]

#output
output1 = open("samplelocimatrix.csv",'w')

#open batch file
with open(file1,'r') as set1:
    catalogdict = {}     #dictionary with: loci_ID List_sample_IDs
    catalogsele = {}     #dictionary for loci shared by minimum of X samples: loci_ID List_sample_IDs
    for i in set1:
        if "#" not in i:
            i = i.split("\t")
            locid = i[2]
            stackids = i[8].split(",")
            stacksidslist = [] #list of samples IDs in the locus

            for j in stackids:
                j = j.split("_")
                isolateid = j[0]
                if isolateid not in stacksidslist:
                    stacksidslist.append(isolateid)

            #feed catalogdict
            stacksidslist.sort(key=int)
            catalogdict[locid] = stacksidslist

            #feed catalogsele
            if len(stacksidslist) >= int(nsamples): #filtering based on the minimum number of samples sharing the locus
                print locid,len(stacksidslist)
                catalogsele[locid] = stacksidslist

set1.close() #close file 1
del catalogdict #not being used

#DICTIONARIES
locipresence = dict() #dictionary with key idstacks, values list of 0 or 1 standing for ausence or presence in loci
sampledata = dict()   #dictionary with key idstacks, values information about sample isolatecode,lineage,isolate,idstacks
locilist = list()     #list of loci

#open file 2
with open(file2,'r') as set2:
    for i in set2:
        i = i.rstrip()
        i = i.split("\t")
        #isolatecode = i[0]
        #lineage = i[1]
        #isolate = i[2]
        idstacks = i[3]
        sampledata[idstacks] = i
        listofpresence = list()
        for j in catalogsele.items(): #verify for each X selected loci if the specific sample (by id stacks) shares the locus
            locus = j[0]
            if locus not in locilist:
                locilist.append(locus)

            representedsamples = j[1]

            if idstacks in representedsamples:
                listofpresence.append(1)
            else:
                listofpresence.append(0)

        locipresence[idstacks] = listofpresence
set2.close() #close file 2
del catalogsele

output1.write("ISOLATE_CODE LINEAGE ISOLATE ID_STACKS SUM_LOCI_LINE "+" ".join(locilist)+"\n");

for i in sampledata.items():
    suma = sum(locipresence[i[0]])
    output1.write(" ".join(i[1])+" "+str(suma)+" "+" ".join(str(x) for x in locipresence[i[0]])+"\n");

locisum = list()
for i in range(len(locilist)):
    samplelocival = []
#    print i
    for j in locipresence.items():
        samplelocival.append(j[1][i])
    locisum.append(sum(samplelocival))
output1.write("    SUM_INDIV_COLUMN "+" ".join(str(x) for x in locisum));
output1.close()

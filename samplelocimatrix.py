#!/usr/bin/python
#Uses batch_1.catalog.tags.tsv output file from Stacks to generate a csv space delimited matrix file where the lines are the samples and columns the loci.
#If sample shares the respective loci the value is 1, whereas for the opposite the value is 0.
#As a second input this script also takes a tab delimited population file <isolatecode><tab><lineage><tab><isolate><tab><idstacks>
#Finally, a minimum number of samples in a locus must be set as parameter

import sys
file1 = sys.argv[1] #batch
file2 = sys.argv[2] #population file <isolatecode><tab><lineage><tab><isolate><tab><idstacks>
nsamples = sys.argv[3] #minimum number of samples in a locus

#output
outputname = "samplelocimatrix.csv"
output1 = open(outputname,'w')

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
    print i

    for j in locipresence.items():
        samplelocival.append(j[1][i])

    locisum.append(sum(samplelocival))

output1.write("    SUM_INDIV_COLUMN "+" ".join(str(x) for x in locisum));

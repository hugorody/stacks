#!/usr/bin/python
#This script is the part 2nd of "From Stacks to PCA" pipeline
#Input: 3 files and 1 string
#------ 1. List of selected Stacks IDs (Population file) with 3 tab separated columns
#------------- 1 column: Individual file name, for example 003
#------------- 2 column: Individual Stacks ID
#------------- 3 column: Individual Code, for example MY_CODE003
#------ 2. List of loci shared among individuals, can be obtained from batch_1.catalog.tags.tsv
#------ 3. Stacks file "batch_1.catalog.tags.tsv"
#------ 4. nameout. a name for the run, for example analise5
#Outputs:
#------ 1. fasta file containing sequences for each individual of concatenated loci from input 2
#------ 2. text file containing the order of concatenated loci

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]
nameout = sys.argv[4]

### OPEN FILES TO CREATE DICTIONARIES
with open(file1,'r') as set1:
    sampledict = dict() #Dict {idstack:'file ID'}
    samplevsindividualids = dict() #Dict {idstack:'sample user code'}
    individuallist = dict()
    for i in set1:
        i = i.split()
        idstack = i[1]
        sampledict[int(idstack)] = i[0]
        samplevsindividualids[int(idstack)] = i[2]
set1.close()

with open(file2,'r') as set2:
    locidict = dict() #Dict {shared loci ID:''}
    for i in set2:
        i = i.rstrip()
        locidict[int(i)] = ''
set2.close()

with open(file3,'r') as set3:
    catalogids = {} #Dict {loci_ID:[List_sample_IDs]}
    catalogloci = {} #Dict {loci_ID:[List_loci_IDs]}
    consensuseq = {} #Dict {loci_ID:consensus_sequence}
    individualloci = {} #Dict {sampleid:[list of loci in the sample that is shared]}

    for i in set3:
        if "#" not in i:
            i = i.split("\t")
            locid = int(i[2]) #catalog shared loci ID

            if locid in locidict: #only if loci is in selected loci input
            #    consensuseq[locid] = i[9] #feeds dictionary with fasta consensus sequence of shared loci
                stackids = i[8].split(",") #takes the column in batch file where sampleID_locusID are separated by comma and split into list
                idssamplelist = [] #list of samples_IDs in the shared loci
                idslocilist = []   #list of locis_IDs in the shared loci

                for j in stackids:
                    j = j.split("_")
                    sampleid = int(j[0])
                    locus_id = int(j[1])

                    if sampleid in sampledict:
                        if sampleid not in individualloci:
                            individualloci[sampleid] = [locus_id]
                        else:
                            addloci = individualloci[sampleid]
                            addloci.append(locus_id)
                            individualloci[sampleid] = addloci

                    #feed lists
                    idssamplelist.append(sampleid)
                    idslocilist.append(locus_id)

                catalogids[locid] = idssamplelist
                catalogloci[locid] = idslocilist
set3.close()

### CREATE CORRELATION AMONG DICTIONARIES AND GENERATE OUTPUTs
individualcatalog = {} #Dict {StackID_LociID: sequence}
for i in individualloci.items(): #{sampleid:[list of loci in the sample that is shared]}
    f1 = open(file3.replace("batch_1.catalog.tags.tsv","")+"joined_"+sampledict[i[0]]+".tags.tsv","r").readlines()
    print "Creating catalog for the file "+str(sampledict[i[0]])+" Stacks ID "+str(i[0])+"."
    for j in f1:
        if "consensus" in j:
            j = j.rstrip()
            j = j.split("\t")
            if int(j[2]) in i[1]:
                individualcatalog[str(i[0])+"_"+j[2]] = j[9]
            #    print str(i[0])+"_"+j[2],j[9]
#print individualcatalog

output1 = open(nameout+"_concater.fasta","w")
for i in individualloci.items(): #{sampleid:[list of loci in the sample that is shared]}
    individualconcatedseq = []
    for j in catalogids.items(): #{loci_ID:[List_sample_IDs]}
        if int(i[0]) in j[1]:
            manytimes = j[1].count(i[0])
            if manytimes == 1: #eliminate duplicated loci
                positionofindividual = j[1].index(i[0])
                sampleid_locid = str(i[0])+"_"+str(catalogloci[j[0]][positionofindividual])
                #print samplevsindividualids[i[0]],"[Stacks ID: "+str(i[0])+", Shared loci ID: "+str(j[0])+", SampleID_LociID:",sampleid_locid+"]"
                individualconcatedseq.append(individualcatalog[sampleid_locid])
            else:
        #        print "duplicated locus"
                individualconcatedseq.append(str(80 * "?"))
        else:
        #    print "missing locus"
            individualconcatedseq.append(str(80 * "?"))

    output1.write(">"+samplevsindividualids[i[0]]+" [Stacks ID "+str(i[0])+"]\n"+"MMM".join(individualconcatedseq)+"\n");
output1.close()

output2 = open(nameout+"_concater_loci_order.fasta","w")
for i in catalogids:
    output2.write(str(i)+"\n");
output2.close()

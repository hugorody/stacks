#!/usr/bin/python

import sys
file1 = sys.argv[1] #ids list  <file ID><tab><stacks/sample ID><tab><sample user code>  note: isolate ID is normally the name of the fastq file
file2 = sys.argv[2] #loci list <catalog_shared_loci>
file3 = sys.argv[3] #batch_1.catalog.tags.tsv

#OPEN FILES TO CREATE DICTIONARIES
with open(file1,'r') as set1:
    sampledict = dict()
    samplevsisolateids = dict()
    for i in set1:
        i = i.split()
        isolate = i[0]
        idstack = int(i[1])
        sampledict[idstack] = isolate
        samplevsisolateids[idstack] = i[2]
set1.close()

with open(file2,'r') as set2:
    locidict = dict()
    for i in set2:
        i = i.rstrip()
        locidict[int(i)] = ''
set2.close()

with open(file3,'r') as set3:
    catalogids = {}     #dictionary with: loci_ID List_sample_IDs
    catalogloci = {}    #dictionary with: loci_ID List_loci_IDs
    consensuseq = {}    #dictionary with: loci_ID consensus_sequence

    for i in set3:
        if "#" not in i:
            i = i.split("\t")
            locid = int(i[2]) #catalog shared loci ID

            if locid in locidict: #only if loci is in selected loci input
                consensuseq[locid] = i[9] #feeds dictionary with fasta consensus sequence of shared loci
                stackids = i[8].split(",") #takes the column in batch file where sampleID_locusID are separated by comma and split into list
                idssamplelist = [] #list of samples_IDs in the shared loci
                idslocilist = []   #list of locis_IDs in the shared loci

                for j in stackids:
                    j = j.split("_")
                    sampleid = int(j[0])
                    locus_id = int(j[1])
                    #feed lists
                    idssamplelist.append(sampleid)
                    idslocilist.append(locus_id)

                catalogids[locid] = idssamplelist
                catalogloci[locid] = idslocilist
set3.close()

#GET DATA
concater = dict() #concatenated fasta file of sample consensus sequences
ordersamples = list() #order of loci being read
for i in locidict: #for each shared_loci in the input list
    output = open(str(i)+".fasta","w") #creates a fasta file
    ordersamples.append(i) #make a record of the order of loci being read
    for j in sampledict: #verifica quais samples possuem aquele loci consultando na lista de samples_IDs do catalogids
        if j in catalogids[i] and int(catalogids[i].count(j)) == 1: #se a sample estiver no loci and no more than one loci per the sample
            position = int(catalogids[i].index(j)) #captura a posicao da sample na lista dos catalogids
            samplelocus = catalogloci[i][position] #utiliza a posicao da sample para pegar sua respectiva locus id no dicionario catalogloci

            f1 = open("joined_"+sampledict[j]+".tags.tsv","r").readlines() #open respective sample.tags.tsv file and read lines
            for line in f1:
                if "consensus" in line: #only take consensus sequences lines
                    line = line.split("\t")
                    if samplelocus == int(line[2]): #if sample locus id is found in line
                        sampleseq = line[9] #take the consensus sequence
                        print ">"+str(j)+"_"+str(samplelocus)+" isolate "+str(samplevsisolateids[j])+" catalog loci "+str(i)+"\n"+sampleseq+"\n"
                        output.write(">"+str(j)+"_"+str(samplelocus)+" isolate "+str(samplevsisolateids[j])+" catalog loci "+str(i)+"\n"+sampleseq+"\n"); #write in the shared loci fasta file, with sampleID_lociID as header

                        if j in concater:  #feeds the concater dict with sampleID as keys and respective consensus sequences as values
                            seqseq = concater[j]+"MMM"+sampleseq
                            concater[j] = seqseq
                        else:
                            concater[j] = sampleseq

        else:
            sampleseq = "?" * 80
            print ">"+str(j)+"_"+str(samplelocus)+" isolate "+str(samplevsisolateids[j])+" catalog loci "+str(i)+"\n"+sampleseq+"\n"
            output.write(">"+str(j)+"_"+str(samplelocus)+" isolate "+str(samplevsisolateids[j])+" catalog loci "+str(i)+"\n"+sampleseq+"\n");

            if j in concater:
                seqseq = concater[j]+"MMM"+sampleseq
                concater[j] = seqseq
            else:
                concater[j] = sampleseq
    output.close()

#concatenated file
output1 = open("concater.fasta","w")
for i in concater.items():
    output1.write(">"+samplevsisolateids[i[0]]+" Stacks_sample_ID "+str(i[0])+"\n"+i[1]+"\n");
output1.close()

#print order of concatenated locis
output2 = open("concater_info.txt","w")
for i in ordersamples:
    output2.write(str(i)+"\n");
output2.close()

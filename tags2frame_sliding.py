#!/usr/bin/python3
#Generate consensus FASTA file from Stacks *.tags.tsv output file
#Blast the FASTA file against CDS reference sequences
#Reads BLAST output file and set new coordinates considering CDS codon frame

import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

input1 = "/home/hugo/Dropbox/Unifesp/Trabalhos/Trabalho1/Blastn2GOtable/3018.tags.tsv"
indiv = "3018"
reference = "./plaza/cds.sbi.csv"
runblast = 0
numberofcopies = 5
mismatches = 10
gapopenings = 10

#GENERATE CONSENSUS FASTA FROM TAGS.TSV
outfasta = open(indiv+".fasta","w")
with open(input1,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "consensus" in i:
            i = i.split("\t")
            locusid = i[1] + "_" + i[2]
            locuseq = i[9]
            outfasta.write(">" + locusid + "\n" + locuseq + "\n")
outfasta.close()

#BLAST STEP
if runblast == 1:
    makeblastdb = subprocess.Popen("makeblastdb -in " + reference + " -dbtype nucl -out " + reference + ".blastdb",shell=True)
    makeblastdb.wait()
    blastn = subprocess.Popen("blastn -query " + indiv + ".fasta -db " + reference + ".blastdb -evalue 1e-20 -outfmt 6 -out " + indiv + ".blastn", shell=True)
    blastn.wait()

#BLAST FILTERING
blastnfilt = {}
with open(indiv + ".blastn","r") as setblastn:
    for i in setblastn:
        i = i.rstrip()
        i = i.split("\t")
        query = i[0]
        refer = i[1] #reference
        identity = float(i[2])
        length = int(i[3])

        if int(i[4]) <= mismatches and int(i[5]) <= gapopenings: #CUTOFF BLAST

            if query not in blastnfilt:
                blastnfilt[query] = i
            else:
                if identity > float(blastnfilt[query][2]) and length > int(blastnfilt[query][3]): #FILTER best hits based on identity and alignment length
                    blastnfilt[query] = i

#PARSE TAGS STACKS FILE
locus = {}
locusblacklist = {}
with open(input1,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "#" not in i and "consensus" not in i and "model" not in i:
            i = i.split("\t")
            locusid = i[1] + "_" + i[2]
            locuseq = i[9]
            readtype = i[6]
            if locusid in blastnfilt and readtype == "primary":
                readindex = int(i[7])
                if locusid not in locus:
                    locus[locusid] = [locuseq]
                    locusblacklist[locusid] = readindex
                else:
                    if int(readindex) > int(locusblacklist[locusid]):
                        addlocuseq = locus[locusid]
                        addlocuseq.append(locuseq)
                        locus[locusid] = addlocuseq
                        locusblacklist[locusid] = readindex

################################################################################
def sliding(startref,stopref):
    if startref % 3 == 0:
        startref = startref
        fatorstart = 0
    else:
        if int(startref + 1) % 3 == 0:
            startref = startref + 1
            fatorstart = 1
        else:
            startref = startref + 2
            fatorstart = 2

    frame = stopref - startref
    frame1 = stopref - 1 - startref
    frame2 = stopref - 2 - startref
    if frame % 3 == 0:
        stopref = stopref
        fatorstop = 0
    else:
        if frame1 % 3 == 0:
            stopref = stopref - 1
            fatorstop = 1
        else:
            stopref = stopref - 2
            fatorstop = 2

    return(startref,stopref,fatorstart,fatorstop)

revcomp = {} #if the gbs sequence is reverse complement 1, or not 0
coordgbs = {} #start stop in gbs
coordref = {} #start stop in ref
dictref = {}
for i in blastnfilt.values():
    query = i[0]
    subject = i[1] #reference
    startquery = int(i[6]) #start of alignment in query
    stopquery = int(i[7]) #end of alignment in query
    startsubject = int(i[8]) #start of alignment in subject/reference
    stopsubject = int(i[9]) #end of alignment in subject/reference

    if len(locus[query]) >= numberofcopies:
        if stopsubject < startsubject: ### MAIN CONDITION
            revcomp[query] = 1 #ALIGNMENT IS REVERSE COMPLEMENT
            startref,stopref,fatorstart,fatorstop = sliding(stopsubject,startsubject)
            startquery = startquery + fatorstop - 1
            stopquery = stopquery - fatorstart - 1

            coordref[subject + "--" + query] = [startref,stopref]
            dictref[subject] = ''
            coordgbs[query] = [startquery,stopquery]

        else: ### NEGATIVE OF MAIN CONDITION
            revcomp[query] = 0 #ALIGNMENT IS NOT REVERSE COMPLEMENT
            startref,stopref,fatorstart,fatorstop = sliding(startsubject,stopsubject)
            startquery = startquery + fatorstart
            stopquery = stopquery - fatorstop

            coordref[subject + "--" + query] = [startref,stopref]
            dictref[subject] = ''
            coordgbs[query] = [startquery,stopquery]

################################################################################
cdsref = {}
for record in SeqIO.parse(reference, "fasta"):
    header = record.id
    if header in dictref:
        cdsref[header] = record.seq

###
for x in blastnfilt.values():
    query = x[0]
    refer = x[1]
    copies = []

    if len(locus[query]) >= numberofcopies: #FILTER BY NUMBER OF COPIES
        output = open(query + "_" + refer + ".gbs.fas","w")
        if revcomp[query] == 1:
            for w in locus[query]:
                seq = w[coordgbs[query][0]:coordgbs[query][1]]
                seq = str(Seq(seq).reverse_complement()) #convert seq to reverse complement
                copies.append(seq)
        else:
            for w in locus[query]:
                seq = w[coordgbs[query][0]:coordgbs[query][1]]
                copies.append(seq)

        #GET SEQUENCE FROM REFERENCE
        myrefseq = cdsref[refer][coordref[refer + "--" + query][0]:coordref[refer + "--" + query][1]]

        #WRITE FASTA OUTPUT FILES
        output.write(">" + refer + " " + str(coordref[refer + "--" + query][0]) + "," + str(coordref[refer + "--" + query][1]) + "\n" + str(myrefseq) + "\n")
        if int(coordref[refer + "--" + query][0]) % 3 != 0: #report if any reference sequence starts out of reading frame
            print ("ERRO in file:",query + "_" + refer + ".gbs.fas","OUT OF FRAME")

        for j in copies:
            output.write(">" + query + " " + str(revcomp[query]) + " " + str(coordgbs[query][0]) + "," + str(coordgbs[query][1]) + "\n" + str(j) + "\n")
        output.close()

#MUSCLE ALIGNMENT STEP
for x in blastnfilt.values():
    query = x[0]
    refer = x[1]
    if len(locus[query]) >= numberofcopies: #FILTER BY NUMBER OF COPIES
        print ("Aligning " + query + "_" + refer)
        muscle = subprocess.Popen("muscle -in " + query + "_" + refer + ".gbs.fas" + " -out " + query + "_" + refer + ".aln.fas",shell=True)
        muscle.wait()

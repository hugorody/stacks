#!/usr/bin/python3
#Generate consensus FASTA file from Stacks *.tags.tsv output file
#Blast the FASTA file against CDS reference sequences
#Reads BLAST output file and set new coordinates considering CDS codon frame

import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

input1 = "3018.tags.tsv"
indiv = "3018"
reference = "cds.sbi.fasta"
runblast = 0 #run 1, not run 0
readlength = 80

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

        if query not in blastnfilt:
            blastnfilt[query] = i
        else:
            if identity > float(blastnfilt[query][2]) and length > int(blastnfilt[query][3]): #filter blastn results getting best hits based on identity and alignment length
                blastnfilt[query] = i

#find reading frame stop alignment in reference
def checkstopref(startref,stopref):
    if (stopref - startref) % 3 == 0:
        stopref = stopref
    else:
        newstop = stopref - 1
        if (newstop - startref) % 3 == 0:
            stopref = newstop
        else:
            stopref = stopref - 2
    return (stopref)

#find reading frame start alignment in reference
def checkstartref(startref):
    newframe = startref + 1
    if newframe % 3 == 0:
        startref = newframe
    else:
        newframe = startref + 2
        if newframe % 3 == 0:
            startref = newframe

    return (startref)

revcomp = {} #if the gbs sequence is reverse complement 1, or not 0
coordgbs = {} #start stop in gbs
coordref = {} #start stop in ref
dictref = {}
for i in blastnfilt.values():
    query = i[0]
    refer = i[1] #reference
    identity = float(i[2])
    length = int(i[3]) #alignment length
    gapopenings = int(5) #number of gap openings
    startquery = int(i[6]) #start of alignment in query
    stopquery = int(i[7]) #end of alignment in query
    startref = int(i[8]) #start of alignment in subject/reference
    stopref = int(i[9]) #end of alignment in subject/reference

    if stopref < startref: ### MAIN CONDITION
        #if so, GBS locus sequence is aligned as reverse complement
        #stopref now represents the start of aligned in reference
        revcomp[query] = 1 #feed dictionary: identity query as reverse complement in dict

        #verify if start of alignment in reference is in reading frame
        frame = stopref - 1
        if frame % 3 == 0: #if the number of nucl in frame is divided into a set of consecutive triplets
            startref = checkstopref(stopref,startref) - 1 #startref works as stopref, uses checkstopref function to find a reading frame stopref
            startquery = readlength - stopquery
            stopquery = startquery + int(startref - stopref) + 1 #query is always aligned in forward. Find stopquery based on the length of reference alignment
            #feed dictionaries
            coordref[refer + "--" + query] = [stopref - 1,startref]
            dictref[refer] = ''
            coordgbs[query] = [startquery,stopquery]

        else: #if the number of nucl in frame is NOT divided into a set of consecutive triplets
            stopref = checkstartref(stopref) #stopref work as startref, uses checkstartref function to find a reading frame startref
            startref = checkstopref(stopref,startref) - 1
            startquery = int(readlength - stopquery) + 1
            stopquery = startquery + int(startref - stopref) + 1
            #feed dictionaries
            coordref[refer + "--" + query] = [stopref,startref + 1] #
            dictref[refer] = ''
            coordgbs[query] = [startquery,stopquery]

    else: ### NEGATIVE OF MAIN CONDITION
        #if so, GBS locus is aligned in forward of reference
        revcomp[query] = 0 #feed dictionary: identity query as NOT reverse complement in dict

        #verify if start of alignment in reference is in reading frame
        frame = startref - 1
        if frame % 3 == 0: #if the number of nucl in frame is divided into a set of consecutive triplets
            startref = startref -1
            stopref = checkstopref(startref,stopref) #uses checkstopref function to find a reading frame stopref
            startquery = startquery - 1
            stopquery = startquery + int(stopref - startref)

            coordref[refer + "--" + query] = [startref,stopref] #
            dictref[refer] = ''
            coordgbs[query] = [startquery,stopquery] #

        else:
            startref = checkstartref(startref) #uses checkstartref function to find a reading frame startref
            stopref = checkstopref(startref,stopref) #uses checkstopref function to find a reading frame stopref
            startquery = startquery + 1
            stopquery = startquery + int(stopref - startref)
            coordref[refer + "--" + query] = [startref,stopref] #
            dictref[refer] = ''
            coordgbs[query] = [startquery,stopquery] #

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
            if locusid in coordgbs and readtype == "primary":
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

cdsref = {}
for record in SeqIO.parse(reference, "fasta"):
    header = record.id
    if header in dictref:
        cdsref[header] = record.seq

###
for x in blastnfilt.values():
    query = x[0]
    refer = x[1] #reference
    copies = []
    #GET SEQUENCES FROM GBS

    if len(locus[query]) >= 5: #FILTER BY NUMBER OF COPIES
        output = open(query + "_" + refer + ".aln.fas","w")
        if int(revcomp[query]) == 1:
            for w in locus[query]:
                seq = str(Seq(w).reverse_complement()) #here has to be reverse complement
                seq = seq[coordgbs[query][0]:coordgbs[query][1]]
                copies.append(seq)
        else:
            for w in locus[query]:
                seq = w[coordgbs[query][0]:coordgbs[query][1]]
                copies.append(seq)

    #GET SEQUENCE FROM REFERENCE
        myrefseq = cdsref[refer][coordref[refer + "--" + query][0]:coordref[refer + "--" + query][1]]

        output.write(">" + refer + " " + str(coordref[refer + "--" + query][0]) + "," + str(coordref[refer + "--" + query][1]) + "\n" + str(myrefseq) + "\n")
        if int(coordref[refer + "--" + query][0]) % 3 != 0: #report if any reference sequence starts out of reading frame
            print ("ERRO in file:",query + "_" + refer + ".aln.fas","OUT OF FRAME")

        for j in copies:
            output.write(">" + query + " " + str(revcomp[query]) + " " + str(coordgbs[query][0]) + "," + str(coordgbs[query][1]) + "\n" + str(j) + "\n")
        output.close()

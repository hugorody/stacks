#!/usr/bin/python3

import sys
import subprocess
from Bio import SeqIO
import os

file1 = sys.argv[1] #concater file
file2 = sys.argv[2] #reference fasta
referencename = sys.argv[3]

#creates seqs dictionary with sequences for each loci
with open(file1,"r") as set1:
    seqs = {}
    for line in set1:
        if line != "\n":
            line = line.rstrip()
            if line[0] != '>':
                line = line.split("MMM")
                count = 0
                for i in line[1:]:
                    count +=1
                    if "?" not in i and str(count) not in seqs:
                        seqs[str(count)] = i
#fasta file with loci sequence that will be used as query during BLAST searches against reference genome
fastaoutput = open("concater_sample.fasta","w")
for i in seqs.items():
    fastaoutput.write(">"+i[0]+"\n"+i[1]+"\n");
fastaoutput.close()

#BLAST steps
makeblastdb = subprocess.Popen("makeblastdb -in " + file2 + " -dbtype nucl -out " + referencename + ".blastdb", shell=True)
makeblastdb.wait()

blastn = subprocess.Popen("blastn -query concater_sample.fasta -db " + referencename + ".blastdb -evalue 1e-20 -outfmt 6 -out " + referencename + ".blastn", shell=True)
blastn.wait()

#read blast output
with open(referencename+".blastn","r") as set2:
    ref_coordinates = {}
    loc_coordinates = {}
    for i in set2:
        i = i.rstrip()
        i = i.split("\t")
        query = i[0]
        ref = i[1]
        startquery = i[6]
        stopquery = i[7]
        coordinatesloc = []
        coordinatesloc.append(startquery)
        coordinatesloc.append(stopquery)
        if query not in loc_coordinates:
            loc_coordinates[query] = coordinatesloc

        startref = i[8]
        stopref = i[9]
        coordinatesref = []
        coordinatesref.append(ref)
        coordinatesref.append(startref)
        coordinatesref.append(stopref)
        if query not in ref_coordinates:
            ref_coordinates[query] = coordinatesref

outputpa = open(referencename+"_loci_presence_absence.txt","w")
for i in seqs:
    if i in ref_coordinates:
        outputpa.write(i+" yes\n");
    else:
        outputpa.write(i+" no\n");
outputpa.close()

#catalog with genome ref sequences
referenceseqs = {}
for record in SeqIO.parse(file2, "fasta"):
    referenceseqs[record.id] = record.seq

finalrefseqs = {}
for i in ref_coordinates.items():
    locus = i[0]
    refheader = i[1][0]
    startref = int(i[1][1])
    stopref = int(i[1][2])

    startloc = int(loc_coordinates[locus][0])
    stoploc = int(loc_coordinates[locus][1])
    variableinstart = (startloc - 1) * "?"
    variableinstop = (80 - stoploc) * "?"
    #print (startloc,stoploc)

    if refheader in referenceseqs:
        if startref < stopref:
            if len(variableinstart+referenceseqs[refheader][startref-1:stopref]+variableinstop) == 80:
                finalrefseqs[locus] = variableinstart+referenceseqs[refheader][startref-1:stopref]+variableinstop
            elif len(variableinstart+referenceseqs[refheader][startref-1:stopref]+variableinstop) < 80:
                addtoend = int(80 - len(variableinstart+referenceseqs[refheader][startref-1:stopref]+variableinstop)) * "?"
                finalrefseqs[locus] = variableinstart+referenceseqs[refheader][startref-1:stopref]+variableinstop+addtoend
        else:
            if len(variableinstart+referenceseqs[refheader][stopref+1:startref]+variableinstop) == 80:
                finalrefseqs[locus] = variableinstart+referenceseqs[refheader][stopref+1:startref].reverse_complement()+variableinstop
            elif len(variableinstart+referenceseqs[refheader][stopref+1:startref]+variableinstop) < 80:
                addtoend = int(80 - len(variableinstart+referenceseqs[refheader][stopref+1:startref]+variableinstop)) * "?"
                finalrefseqs[locus] = variableinstart+referenceseqs[refheader][stopref+1:startref].reverse_complement()+variableinstop+addtoend

#write to concatenated file given as input the reference sequences of loci concatenated
concatenatedfile = open(file1,"a")
concatenatedfile.write(">"+referencename+"\n");
presence_absense_inref = {} #dict with key locus number and presence absence information as values
for i in range(1,len(seqs)+1):
    if str(i) in finalrefseqs:
        presence_absense_inref[i] = "presence"
        if len(finalrefseqs[str(i)]) < 80:
            print (i,len(finalrefseqs[str(i)]))
        concatenatedfile.write("MMM"+str(finalrefseqs[str(i)]));
    else:
        concatenatedfile.write("MMM" + 80 * "?");
        presence_absense_inref[i] = "absence"
concatenatedfile.write("\n");
concatenatedfile.close()

#read concater file again but now it has the genome reference sequences
with open(file1,"r") as set1:
    seqsconcat = {}
    for line in set1:
        if line != "\n":
            line = line.rstrip()
            if line[0] == '>':

                name=line[1:]
                seqsconcat[name] = []
            else:
                line = line.split("MMM")
                line = line[1:]
                seqsconcat[name] = line


outputfinal = open("concater_loci_plus_ref.fasta","w")
for i in seqsconcat.items():
    outputfinal.write(">"+i[0]+"\n");
    for j in presence_absense_inref.items():
        if j[1] == "presence":
            position = int(j[0]) - 1
            outputfinal.write("MMM"+str(i[1][position]));
    outputfinal.write("\n");
outputfinal.close()

#remove temporary files
os.remove(referencename+".blastdb.nhr")
os.remove(referencename+".blastdb.nin")
os.remove(referencename+".blastdb.nsq")
#os.remove("concater_sample.fasta")

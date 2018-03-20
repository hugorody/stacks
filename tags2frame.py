#!/usr/bin/python3
#Generate consensus FASTA file from Stacks *.tags.tsv output file
#Blast the FASTA file against CDS reference sequences
#Reads BLAST output file and set new coordinates considering CDS codon frame

import subprocess

input1 = "/home/hugo/Dropbox/Unifesp/Trabalhos/Trabalho1/Blastn2GOtable/3018.tags.tsv"
indiv = "3018"
reference = "./plaza/cds.sbi.csv"
runblast = 0

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

revcomp = {} #if the gbs sequence is reverse complement 1, or not 0
coordgbs = {} #start stop in gbs
coordref = {} #start stop in ref
for i in blastnfilt.values():
    query = i[0]
    refer = i[1] #reference
    identity = float(i[2])
    length = int(i[3])
    gapopenings = int(5)
    startquery = int(i[6])
    stopquery = int(i[7])
    startref = int(i[8])
    stopref = int(i[9])
    if stopref < startref: #pegar reversa complementar do gbs
        #stopref passara a ser o inicio do alinhamento do gbs reverso complementar na referencia
        revcomp[query] = 1 #adiciona como revcomp ao dicionario
        frame = stopref - 1 #verifica o n de nucl antes do inicio do alinhamento
        if frame % 3 == 0: #se o inicio do alinhamento ja esta em frame em relacao a referencia
            startref = checkstopref(stopref,startref)
            coordref[refer] = [stopref,startref] #FRAME VERIFICADO primeiro valor relativo ao inicio do alinhamento

            stopquery = startquery + int(startref - stopref)
            coordgbs[query] = [startquery,stopquery]


        else: #se inicio do alinhamento nao esta em frame em relacao a referencia
            if stopref % 3 == 0:
                stopref = stopref + 1
            else:
                newframe = stopref + 1
                if newframe % 3 == 0:
                    stopref = newframe
                else:
                    stopref = stopref + 2

            startref = checkstopref(stopref,startref)
            coordref[refer] = [stopref,startref] #FRAME VERIFICADO

            stopquery = startquery + int(startref - stopref)
            coordgbs[query] = [startquery,stopquery]


    else: #pega a sequencia do gbs normal
        revcomp[query] = 0 #adiciona como nao revcomp ao dicionario
        frame = startref - 1
        if frame % 3 == 0: #se o inicio do alinhamento ja esta em frame em relacao a referencia
            stopref = checkstopref(startref,stopref) #use function to gather stopref
            coordref[refer] = [startref,stopref] #FRAME VERIFICADO

            stopquery = startquery + int(stopref - startref)
            coordgbs[query] = [startquery,stopquery] #FRAME VERIFICADO

        else:
            if startref % 3 ==0:
                startref = startref + 1
            else:
                newframe = startref + 1
                if newframe % 3 == 0:
                    startref = newframe
                else:
                    startref = startref + 2

            stopref = checkstopref(startref,stopref) #use function to gather stopref
            coordref[refer] = [startref,stopref] #FRAME VERIFICADO

            stopquery = startquery + int(stopref - startref)
            coordgbs[query] = [startquery,stopquery] #FRAME VERIFICADO

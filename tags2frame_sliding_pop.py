#!/usr/bin/python3

from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio import SeqIO
import multiprocessing
from multiprocessing import Process, Manager
import subprocess
import itertools
import glob,os
import time

print ("Start of run:",time.ctime())

#variables
fileparent1 = "3018"
fileparent2 = "3046"
batch = "batch_1.catalog.tags.tsv"
workdir = "/storage2/hugo/sugarcane/denovo_out09/" #directory of Stacks tags.tsv files
reference = "/storage1/referencegenomes/cds.sbi.csv"
copiescutoff = 2 #minimum number of Stacks for a Locus Stack in progenitors
num_cores = 30
runblast = 1 #1 for yes, 0 for no
mismatches = 10
gapopenings = 10


outlog = open(workdir + "tags2frame.log","w")
################################################################################
#Relation between files and Stacks IDs
print ("Creating relationship between tags.tsv files and its correspondents Stacks IDs")
outlog.write("Creating relationship between tags.tsv files and its correspondents Stacks IDs "+str(time.ctime())+"\n")
def getstackID(filetagstsv):
    with open(filetagstsv,"r") as settagstsv:
        for i in itertools.islice(settagstsv, 1, 2):
            i = i.rstrip()
            i = i.split("\t")
            stackID = int(i[1])
    return(stackID)

filestackIDs = {}
os.chdir(workdir) #list all files in the input directory 1
#print ("Files tags.tsv in working directory:")
for file in glob.glob("*.tags.tsv"): #and filter by *.fq
    if file != batch:
        namefile = str(file).replace(".tags.tsv","")

        if fileparent1 in file:
            parent1 = getstackID(file) #DEFINE VARIABLE parent1
            filestackIDs[parent1] = namefile
            dictname = "dict" + namefile
            dictname = {}
        elif fileparent2 in file:
            parent2 = getstackID(file) #DEFINE VARIABLE parent2
            filestackIDs[parent2] = namefile
        else:
            idfile = getstackID(file)
            filestackIDs[idfile] = namefile

################################################################################
print ("Defining number of Stacks Locus copies in progenitors.")
outlog.write("Defining number of Stacks Locus copies in progenitors "+str(time.ctime())+"\n")
#Dicts number of COPIES of progenitors
ncopiesPARENT1 = {} #{stackID_locusID:int(ncopies)}
ncopiesPARENT2 = {} #{stackID_locusID:int(ncopies)}

#PARSE both progenitors STACKS files
with open(filestackIDs[parent1]+".tags.tsv","r") as set1:
    for i in set1:
        i = i.rstrip()
        if "#" not in i and "consensus" not in i and "model" not in i:
            i = i.split("\t")
            locusid = int(i[2])
            readtype = i[6]
            if readtype == "primary":
                readindex = int(i[7]) + 1
                ncopiesPARENT1[locusid] = readindex

with open(str(filestackIDs[parent2])+".tags.tsv","r") as set2:
    for i in set2:
        i = i.rstrip()
        if "#" not in i and "consensus" not in i and "model" not in i:
            i = i.split("\t")
            locusid = int(i[2])
            readtype = i[6]
            if readtype == "primary":
                readindex = int(i[7]) + 1
                ncopiesPARENT2[locusid] = readindex

################################################################################
print ("Creating dictionaries and FASTA file for catalog.")
outlog.write("Creating dictionaries and FASTA file for catalog "+str(time.ctime())+"\n")
#Dicts of loci shared by progenitors and individuals
sharedPARENT12 = {} #{catalogID:stacksID1_locusID,stacksID2_locusID,...} loci shared by 3018 and 3046
sharedPARENT1 = {} #{catalogID:stacksID1_locusID,stacksID2_locusID,...} loci not present in 3046
sharedPARENT2 = {} #{catalogID:stacksID1_locusID,stacksID2_locusID,...} loci not present in 3018
indivloci = {} #{stackindivID:[stacklocusID1,stacklocusID2,...]}
outfasta = open("catalogseq.fasta","w")
with open(workdir+batch,"r") as setbatch:
    for i in setbatch:
        i = i.rstrip()
        if "#" not in i:
            i = i.split("\t")
            catalogID = i[2]
            catalogSHARING = i[8].split(",") #list of stacksID_locusID
            catalogSEQ = i[9]

            #HERE CAN BE SET A MINIMUM LEN FOR catalogSHARING AS A FILTER

            liststacksIDs = []
            listlocusIDs = []
            for j in catalogSHARING:
                j = j.split("_")
                liststacksIDs.append(int(j[0]))
                listlocusIDs.append(j[1])

            if int(parent1) in liststacksIDs and int(parent2) in liststacksIDs: #IF LOCUS SHARED BY PARENTs 1 AND 2
                locusparent1 = int(listlocusIDs[liststacksIDs.index(int(parent1))])
                locusparent2 = int(listlocusIDs[liststacksIDs.index(int(parent2))])
                if ncopiesPARENT1[locusparent1] >= copiescutoff and ncopiesPARENT2[locusparent2] >= copiescutoff:
                    sharedPARENT12[catalogID] = catalogSHARING
                    outfasta.write(">" + catalogID + "\n" + catalogSEQ + "\n")

                    for x in catalogSHARING: #FEED DICT INDIVLOCI
                        x = x.split("_")
                        if int(x[0]) not in indivloci:
                            indivloci[int(x[0])] = [int(x[1])]
                        else:
                            addlocus = indivloci[int(x[0])]
                            addlocus.append(int(x[1]))
                            indivloci[int(x[0])] = addlocus

            elif int(parent1) in liststacksIDs and int(parent2) not in liststacksIDs: #IF LOCUS ONLY IN PARENT 1
                locusparent1 = int(listlocusIDs[liststacksIDs.index(int(parent1))])
                if ncopiesPARENT1[locusparent1] >= copiescutoff:
                    sharedPARENT1[catalogID] = catalogSHARING
                    outfasta.write(">" + catalogID + "\n" + catalogSEQ + "\n")

                    for x in catalogSHARING: #FEED DICT INDIVLOCI
                        x = x.split("_")
                        if int(x[0]) not in indivloci:
                            indivloci[int(x[0])] = [int(x[1])]
                        else:
                            addlocus = indivloci[int(x[0])]
                            addlocus.append(int(x[1]))
                            indivloci[int(x[0])] = addlocus

            elif int(parent1) not in liststacksIDs and int(parent2) in liststacksIDs: #IF LOCUS ONLY IN PARENT 2
                locusparent2 = int(listlocusIDs[liststacksIDs.index(int(parent2))])
                if ncopiesPARENT2[locusparent2] >= copiescutoff:
                    sharedPARENT2[catalogID] = catalogSHARING
                    outfasta.write(">" + catalogID + "\n" + catalogSEQ + "\n")

                    for x in catalogSHARING: #FEED DICT INDIVLOCI
                        x = x.split("_")
                        if int(x[0]) not in indivloci:
                            indivloci[int(x[0])] = [int(x[1])]
                        else:
                            addlocus = indivloci[int(x[0])]
                            addlocus.append(int(x[1]))
                            indivloci[int(x[0])] = addlocus

outfasta.close()
print ("Number of individuals in population:",len(indivloci))
outlog.write("Number of individuals in population: "+str(len(indivloci))+", "+str(time.ctime())+"\n")
################################################################################
print ("Starting BLAST step.")
outlog.write("Starting BLAST step "+str(time.ctime())+"\n")

# BLAST AGAINST REFERENCE

if runblast == 1:
    makeblastdb = subprocess.Popen("makeblastdb -in " + reference + " -dbtype nucl -out " + reference + ".blastdb",shell=True)
    makeblastdb.wait()
    blastn = subprocess.Popen("blastn -query catalogseq.fasta -db " + reference + ".blastdb -evalue 1e-20 -outfmt 6 -out catalogseq.blastn", shell=True)
    blastn.wait()

#FILTERED RESULTS FROM BLASTN
blastnfilt = {} #FILTER BY QUERY
with open("catalogseq.blastn","r") as setblastn:
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


################################################################################
#DICT CONTAINING GBS LOCUS COPIES SEQUENCES

def createcatalog(stackindivID,stackslocusIDs): #from dict indivloci {stackindivID:[stacklocusID1,stacklocusID2,...]}
    locusblacklist = {}
    passing = 0
    if stackindivID in filestackIDs: #RUN FUNCTION ONLY FOR FILES IN workdir
        print ("Creating catalog for:",stackindivID,filestackIDs[stackindivID])
        outlog.write("Creating catalog for: "+str(stackindivID)+" "+str(filestackIDs[stackindivID])+" "+str(time.ctime())+"\n")
        with open(workdir+str(filestackIDs[stackindivID])+".tags.tsv","r") as set1: #open file stackindivID
            for j in set1:
                j = j.rstrip()
                j = j.split("\t")
                if "#" not in j:
                    if "consensus" in j:
                        if int(j[2]) in stackslocusIDs:
                            passing = 1
                    #j = j.split("\t")
                    #locusid = int(j[2])
                    #readtype = j[6]
                    #print (passing,j)
                    if passing == 1 and "primary" in j:
                        readindex = int(j[7])
                        #print (passing,readindex,j)
                        #locuseq = j[9]
                        stacksID_locusID = str(stackindivID) + "_" + str(j[2])
                        if stacksID_locusID not in stackIDlocicopies:
                            stackIDlocicopies[stacksID_locusID] = [j[9]]
                            locusblacklist[stacksID_locusID] = readindex
                            #print (stacksID_locusID,j[9])
                        else:
                            if readindex > locusblacklist[stacksID_locusID]:
                                addlocuseq = stackIDlocicopies[stacksID_locusID]
                                addlocuseq.append(j[9])
                                stackIDlocicopies[stacksID_locusID] = addlocuseq
                                locusblacklist[stacksID_locusID] = readindex
                                #print (stacksID_locusID,addlocuseq)
                    if "secondary" in j:
                        passing = 0

if os.path.isfile(workdir + "stackIDlocicopiesoutput.tmp") is True: #if the file already exist
    print ("Temporary file found. Parsing to create stackIDlocicopies dict.")
    stackIDlocicopies = {}
    with open(workdir + "stackIDlocicopiesoutput.tmp","r") as sett1:
        for i in sett1:
            i = i.rstrip()
            i = i.split(" ")
            stackIDlocicopies[i[0]] = []
            listseq = i[1].split(",")
            for j in listseq:
                addj = stackIDlocicopies[i[0]]
                addj.append(j)
                stackIDlocicopies[i[0]] = addj
else:
    print ("Temporary file NOT FOUND.")
    manager = Manager()
    stackIDlocicopies = manager.dict()
    Parallel( n_jobs=num_cores) ( delayed(createcatalog) (i[0],i[1]) for i in indivloci.items() )

    #CREATE Temporary FILE
    stackIDlocicopiesoutput = open("stackIDlocicopiesoutput.tmp","w")
    for i in stackIDlocicopies.items():
        stackIDlocicopiesoutput.write(str(i[0]) + " " + ",".join(i[1]) + "\n")
    stackIDlocicopiesoutput.close()

################################################################################
#SLIDING FRAME
print ("Slinding Frame step")
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

revcomp = {} #{catalogID:1} if the gbs sequence is reverse complement 1, or not 0
coordgbs = {} #{catalogID:[int(start),int(stop)]} start stop in gbs
coordref = {} #{subject:[int(start),int(stop)]} start stop in ref
dictref = {}

for i in blastnfilt.values():
    query = i[0]
    subject = i[1] #reference
    startquery = int(i[6]) #start of alignment in query
    stopquery = int(i[7]) #end of alignment in query
    startsubject = int(i[8]) #start of alignment in subject/reference
    stopsubject = int(i[9]) #end of alignment in subject/reference

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
#PARSE REFERENCE FASTA
print ("Parsing Reference FASTA file.")
cdsref = {}
for record in SeqIO.parse(reference, "fasta"):
    header = record.id
    if header in dictref:
        cdsref[header] = record.seq
################################################################################
#FINAL STEP
print ("Final step")
for x in blastnfilt.values():
    query = x[0]
    refer = x[1]

    if query in sharedPARENT12: #filter
        headers = []
        copies = []
        output = open(query + "_" + refer + ".gbs.fas","w")

        if revcomp[query] == 1:
            for w in sharedPARENT12[query]:
                for stackcopy in stackIDlocicopies[w]:
                    seq = stackcopy[coordgbs[query][0]:coordgbs[query][1]]
                    seq = str(Seq(seq).reverse_complement()) #convert seq to reverse complement
                    headers.append(w)
                    copies.append(seq)
        else:
            for w in sharedPARENT12[query]:
                for stackcopy in stackIDlocicopies[w]:
                    seq = stackcopy[coordgbs[query][0]:coordgbs[query][1]]
                    headers.append(w)
                    copies.append(seq)

        #GET SEQUENCE FROM REFERENCE
        myrefseq = cdsref[refer][coordref[refer + "--" + query][0]:coordref[refer + "--" + query][1]]

        #WRITE FASTA OUTPUT FILES
        output.write(">" + refer + "\n" + str(myrefseq) + "\n")
        if int(coordref[refer + "--" + query][0]) % 3 != 0: #report if any reference sequence starts out of reading frame
            print ("ERRO in file:",query + "_" + refer + ".gbs.fas","OUT OF FRAME")

        count = 0
        for j in copies:
            output.write(">" + headers[count] + "_seq" + str(count) + "_loc" + query + "\n" + str(j) + "\n")
            count += 1

        output.close()

outlog.close()
print ("End of run:",time.ctime())

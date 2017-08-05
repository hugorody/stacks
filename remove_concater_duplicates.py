#!/usr/bin/python
#read concater fasta and remove loci where duplicates has been identified

import sys
file1 = sys.argv[1] #concater fasta
file2 = sys.argv[2] #concater order
file3 = sys.argv[3] #concater duplicates

with open(file1,'r') as f:
    seqs={}
    for line in f:
        line = line.rstrip()
        if line[0] == '>':
            name=line
            seqs[name]=''
        else:
            seqs[name] = seqs[name] + line

    for i in seqs.items():
        indiv = i[0]
        sequences = i[1].split("MMM")
        nloci = len(sequences)
        seqs[indiv] = sequences

with open(file2,'r') as set2:
    order = 0
    lociorder = {}
    for i in set2:
        i = i.rstrip()
        order += 1
        lociorder[i] = order

with open(file3,'r') as set3:
    duplicates = {}
    for i in set3:
        i = i.split(" ")
        duplicates[i[0]] = ''

for i in duplicates:
    del lociorder[i]

#output
output1name = file1.replace(".fasta","_noduplicates.fasta")
output2name = file1.replace(".fasta","_noduplicates_lociorder.txt")

output2 = open(output2name,"w")
selectedloci = {}
for j in seqs.items():
    selectedloci[j[0]] = ''
    for i in lociorder.items():
        output2.write(str(i[0])+"\n");
        position = i[1] - 1
        selectedloci[j[0]] = selectedloci[j[0]] + "MMM" + j[1][position]
output2.close()

output1 = open(output1name,"w")
for i in selectedloci.items():
    output1.write(i[0]+"\n"+i[1]+"\n");
output1.close()

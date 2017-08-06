#!/usr/bin/python
#read concater fasta and remove loci where duplicates has been identified
#generates two outputs: 1 is the space separated order of individuals and the 2 is the dna sequence consecutive in each line

import sys
file1 = sys.argv[1] #concater fasta

with open(file1,'r') as f:
    seqs={}
    for line in f:
        line = line.rstrip()
        if line[0] == '>':
            words=line.split()
            name=words[0][1:]
            seqs[name]=''
        else:
            seqs[name] = seqs[name] + line.replace("MMM","")

transposed = {}
transposedheader = {}
for i in seqs.items():
    header = i[0]
    mylist = list(i[1])
    sizemylist = len(mylist)

    if "header" not in transposedheader:
        transposedheader["header"] = [header]
    else:
        newheader = transposedheader["header"]
        newheader.append(header)
        transposedheader["header"] = newheader


    for j in range(sizemylist):
        if j not in transposed:
            transposed[j] = [mylist[j]]
        else:
            newletter = transposed[j]
            newletter.append(mylist[j])
            transposed[j] = newletter

output1name = file1.replace(".fasta","_isolate_order.csv")
output1 = open(output1name,"w")
for i in transposedheader.values():
    output1.write(" ".join(i));
output1.close()


output2name = file1.replace(".fasta","_iNPUTE.csv")
output2 = open(output2name,"w")
for i in transposed.values():
    output2.write(" ".join(i)+"\n");
output2.close()

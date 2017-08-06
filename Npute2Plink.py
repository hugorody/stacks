#!/usr/bin/python

import sys
file1 = sys.argv[1] #NPUTE imputation output
file2 = sys.argv[2] #order of individuals space delimited generated by Concater2Npute.py

with open(file1,'r') as set1:
    count = 0
    seqs = {}
    for i in set1:
        i = i.rstrip()
        i = i.split(",")
        count += 1

        if "" not in i: #if imputation has not failed
            for j in range(len(i)): #for the number of columns
                if j not in seqs:
                    seqs[j] = [i[j]]
                else:
                    newseq = seqs[j]
                    newseq.append(i[j])
                    seqs[j] = newseq
#        else:
#            print count

with open(file2,'r') as set2:
    for i in set2:
        i = i.split(" ")
        isolates = i

for i in seqs.items():
    print isolates[i[0]],"1"," ".join(i[1]).replace("A","1").replace("C","2").replace("G","3").replace("T","4")
    #print ">"+isolates[i[0]]+"\n"+"".join(i[1])

#!/usr/bin/python
#usage: python script.py file.txt

import matplotlib.pyplot as plt
from pylab import *
import sys

filetxt = sys.argv[1]
outputdir = sys.argv[2]

#open files
try:
    ftxt = open(filetxt)
except IOError:
    print ("Files doesn't exist!")

reads = [] #number of reads is x axis
loci = [] #number of loci is y axis

for i in ftxt:
    i = i.split(" ")
    reads.insert(0, float(i[1]))
    loci.insert(0, float(i[2]))

#calc of line fit
(m,b)=polyfit(reads,loci,1)
yp = polyval([m,b],reads)

plt.plot(reads,yp) #plot line fit

plt.scatter(reads,loci) #plot scatter

plt.xlim(0,)
plt.ylim(0,)
plt.title('Correlation among individuals n reads and predicted n loci')
plt.xlabel('Number of reads')
plt.ylabel('Number of loci')
plt.savefig(outputdir+"nreadsVSnloci_scatterplot.png") #save file
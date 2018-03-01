#!/usr/bin/python3

import pandas
import pylab as pl
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

#READ TAGS.TSV AND GENERATE INPUT FOR KMEANS
with open("3018.tags.tsv","r") as set1:
    indivCATALOG = {} #Dict {indivID:'locusID',}
    locusCATALOG = {} #Dict {locusID:['sequence1','sequence2']}
    locusCATALOGfilter = {} #Dict {locusID:number}
    for i in set1:
        if "#" not in i and "consensus" not in i and "secondary" not in i and "model" not in i:
            i = i.rstrip()
            i = i.split("\t")
            indivID = i[1]
            locusID = i[2]
            locusQUAL = i[6] #primary or secondary
            locusQUALn = i[7] #primary 1 or 2 or 3 ...
            locusSEQ = i[9] #nucl sequence

            if locusID not in locusCATALOG:
                locusCATALOG[locusID] = [locusSEQ]
                locusCATALOGfilter[locusID] = locusQUALn
            else:
                if int(locusQUALn) > int(locusCATALOGfilter[locusID]):
                    addsequence = locusCATALOG[locusID]
                    addsequence.append(locusSEQ)
                    locusCATALOG[locusID] = addsequence
                    locusCATALOGfilter[locusID] = locusQUALn

kmeansinput = open("kmeansinput.tsv","w")
kmeansinput.write("locus,copies\n");

for i in locusCATALOG.items():
    kmeansinput.write(str(i[0]) + "," + str(len(i[1])) + "\n");


################################################################################

variables = pandas.read_csv('kmeansinput.tsv')
X = variables[['copies']]
Y = variables[['locus']]


#calculate k by elbow method
Nc = range(1, 20)
kmeans = [KMeans(n_clusters=i) for i in Nc]
#score = [kmeans[i].fit(X).score(X) for i in range(len(kmeans))]
score = []
clusterscore = {} #Dict {Cluster:Score}
for i in range(len(kmeans)):
    score.append(kmeans[i].fit(X).score(X))
    clusterscore[i] = kmeans[i].fit(X).score(X)

pl.plot(Nc,score)
pl.xlabel('Number of Clusters')
pl.ylabel('Score')
pl.title('Elbow Curve')
pl.show()

#KMEANS
pca = PCA(n_components=1).fit(X)
pca_d = pca.transform(X)
pca_c = pca.transform(Y)

kmeans=KMeans(n_clusters=5) #number of clusters
kmeansoutput=kmeans.fit(X)
kmeansoutput

pl.figure('5 Cluster K-Means')
pl.scatter(pca_c[:, 0], pca_d[:, 0], c=kmeansoutput.labels_)
pl.xlabel('Locus')
pl.ylabel('Number of copies')
pl.title('5 Cluster K-Means')
pl.show()

#output
clusteroutput = open("cluster_output.csv","w")
clusteroutput.write("cluster copies locus\n");
clusterranges = {} #Dict {Cluster:[1,2,3,4]}
clusterloci = {} #Dict {Cluster:[locus1,locus2]}
count = 0
for i in kmeans.predict(X):
    cluster = str(i)
    copies = str(X.values[count]).replace("[","").replace("]","")
    locus = str(Y.values[count]).replace("[","").replace("]","")
    clusteroutput.write(cluster + " " + copies + " " + locus + "\n");

    if cluster not in clusterloci:
        clusterloci[cluster] = [int(locus)]
    else:
        addlocus = clusterloci[cluster]
        addlocus.append(int(locus))
        clusterloci[cluster] = addlocus


    if cluster not in clusterranges:
        clusterranges[cluster] = [int(copies)]
    else:
        addval = clusterranges[cluster]
        addval.append(int(copies))
        clusterranges[cluster] = addval

    count += 1

for i in clusterranges.items():
    valmax = max(i[1])
    valmin = min(i[1])
    print ("Cluster",i[0],"range:",valmin,"to",valmax,"(",len(i[1]),"copies )")

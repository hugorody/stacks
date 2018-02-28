#!/usr/bin/python3

#Utiliza comparações BLAST dos locus preditos pelo Stacks contra as sequências
#de CDS de um genoma do Plaza 4.0. Para cada locus do arquivo output do BLAST,
#são atribuídas as N anotações GO respectivas às X sequências do genoma
#que passaram o filtro do BLAST.

#A primeira coluna da tabela (delimitada por vírgula) corresponde aos locus
#preditos, e as demais correspondem a cada uma das categorias GO com 1 para
#presença e 0 para ausência.

outputname = "locus_go_Stacksdenovo.csv"

with open("zma.blastn","r") as set1:
    locusMATCH = {} #Dict {locus:[match1,mathc2]}
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        if i[0] not in locusMATCH:
            locusMATCH[i[0]] = [i[1]]
        else:
            addlocusmatch = locusMATCH[i[0]]
            addlocusmatch.append(i[1])
            locusMATCH[i[0]] = addlocusmatch

with open("go.zma.csv","r") as set2:
    matchGO = {} #Dict {match:[go1,go2]}
    for i in set2:
        if "#" not in i:
            i = i.rstrip()
            i = i.replace("\"","")
            i = i.split(";")
            if i[0] not in matchGO:
                matchGO[i[0]] = [i[2]]
            else:
                addGO = matchGO[i[0]]
                addGO.append(i[2])
                matchGO[i[0]] = addGO

listofGO = []

for i in locusMATCH.items():
    for j in i[1]:
        if j in matchGO:
            for n in matchGO[j]:
                if n not in listofGO:
                    listofGO.append(n)

locusGO = {} #Dict {locus:[go1,go2...]}
for i in locusMATCH.items():
    for j in i[1]:
        if j in matchGO:
            if i[0] not in locusGO:
                locusGO[i[0]] = matchGO[j]
            else:
                addlocusGO = locusGO[i[0]]
                addlocusGO.append(matchGO[j])
                locusGO[i[0]] = addlocusGO


locusGO_0and1 = {}
for i in locusGO.items():
    listOFlocusGOs = i[1]
    list0and1 = []
    for j in listofGO:
        if j in listOFlocusGOs:
            list0and1.append("1")
        else:
            list0and1.append("0")
    locusGO_0and1[i[0]] = list0and1

output1 = open(outputname,"w")
output1.write("locus," + ",".join(listofGO) + "\n");
for i in locusGO_0and1.items():
    output1.write(str(i[0]) + "," + ",".join(i[1]) + "\n");
output1.close()

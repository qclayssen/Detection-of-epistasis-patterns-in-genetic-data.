#!/usr/bin/env python3
# coding: utf-8

import re
import os
import random
import math
import sys
import warnings

outFolder=str(sys.argv[1])
nbFiles=int(sys.argv[2])
prefix=str(sys.argv[3])
nbVar=int(sys.argv[4])
nbPositive=int(sys.argv[5])
nbNegative=int(sys.argv[6])
nbTotal=nbPositive+nbNegative

try:
    os.mkdir(outFolder)
    print("Dossier de sortie" , outFolder ,  " créé")
except FileExistsError:
    print("Le dossier de sortie" , outFolder ,  " existe déjà")

random.seed()


phiList=[]
globalList=[]
phenoPosList=[]
phenoNegList=[]

while (phenoPosList==[] or phenoNegList==[]):
    alpha=random.uniform(-1,1)
    beta1=random.uniform(-1,1)
    beta2=random.uniform(-1,1)
    beta12=random.uniform(-1,1)


    for X1 in range (0,3):
    	for X2 in range (0,3):
    		value=alpha+beta1*X1+beta2*X2+beta12*X1*X2
    		phiList.append(value)
    		theMiddle=[X1,X2,value]
    		globalList.append(theMiddle)


    minim=min(phiList)
    maxim=max(phiList)


    for i in range(0,len(phiList)):
    	globalList[i][2]=1/(1+(math.exp(-phiList[i])))
    	print(globalList[i][2])
    	if (globalList[i][2]>0.5) :
    		phenoPosList.append([globalList[i][0],globalList[i][1]])
    	else :
    		phenoNegList.append([globalList[i][0],globalList[i][1]])

for i in range(1,nbFiles+1):
    newFileGeno=open("{0}/{1}_Genotype_{2}.csv".format(outFolder,prefix,i),"w")
    newFilePheno=open("{0}/{1}_Phenotype_{2}.csv".format(outFolder,prefix,i),"w")
    newFilePheno.write("Class")
    header=""
    for j in range(1,nbVar-1):
        header=header+"N"+str(j)+","
    header=header+"CAUS1,CAUS2\n"
    line=""
    countPosNeg=0
    for k in range(nbTotal):
        for j in range(nbVar-2):
            line=line+str(random.randrange(3))+","
        if (countPosNeg<nbPositive):
            whichDuo=random.randrange(len(phenoPosList))
            line=line+str(phenoPosList[whichDuo][0])+","+str(phenoPosList[whichDuo][1])
            newFilePheno.write("1\n")
        else :
            whichDuo=random.randrange(len(phenoNegList))
            line=line+str(phenoNegList[whichDuo][0])+","+str(phenoNegList[whichDuo][1])
            newFilePheno.write("0\n")
        line=line+"\n"
        countPosNeg+=1
    newFileGeno.write(header)
    newFileGeno.write(line)
    newFileGeno.close
    newFilePheno.close

print("\n",alpha)
print(beta1)
print(beta2)
print(beta12)
print("Cas",phenoPosList)
print("Témoins",phenoNegList)

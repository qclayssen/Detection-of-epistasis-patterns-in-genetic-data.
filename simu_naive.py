#!/usr/bin/env python3
# coding: utf-8

#//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
#//Naive simulation of Epistasis
#//Created :09/11/18
#//Modified :11/02/2019

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
patternSize=int(sys.argv[7])
scramblePercent=int(sys.argv[8])
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

if (patternSize==2):
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

        for i in range(0,len(phiList)):
        	globalList[i][2]=1/(1+(math.exp(-phiList[i])))
        	if (globalList[i][2]>0.5) :
        		phenoPosList.append([globalList[i][0],globalList[i][1]])
        	else :
        		phenoNegList.append([globalList[i][0],globalList[i][1]])

    for i in range(1,nbFiles+1):
        newFileGeno=open("{0}/{1}_Genotype_{2}.csv".format(outFolder,prefix,i),"w")
        newFilePheno=open("{0}/{1}_Phenotype_{2}.csv".format(outFolder,prefix,i),"w")
        newFilePheno.write("Class\n")
        header=""
        for j in range(1,nbVar-1):
            header=header+"N"+str(j)+","
        header=header+"M0P1,M0P2\n"
        line=""
        countPosNeg=0
        incr=0
        for k in range(nbTotal):
            for j in range(nbVar-2):
                line=line+str(random.randrange(3))+","
            if (countPosNeg<nbPositive):
                scramble=random.randrange(100)
                if (scramble<(scramblePercent)):
                    whichDuo=random.randrange(len(phenoNegList))
                    line=line+str(phenoNegList[whichDuo][0])+","+str(phenoNegList[whichDuo][1])
                    incr+=1
                else :
                    whichDuo=random.randrange(len(phenoPosList))
                    line=line+str(phenoPosList[whichDuo][0])+","+str(phenoPosList[whichDuo][1])
                newFilePheno.write("1\n")
            else :
                if (incr>0):
                    whichDuo=random.randrange(len(phenoPosList))
                    line=line+str(phenoPosList[whichDuo][0])+","+str(phenoPosList[whichDuo][1])
                    incr-=1
                else :
                    whichDuo=random.randrange(len(phenoNegList))
                    line=line+str(phenoNegList[whichDuo][0])+","+str(phenoNegList[whichDuo][1])
                newFilePheno.write("0\n")
            line=line+"\n"
            countPosNeg+=1
        newFileGeno.write(header)
        newFileGeno.write(line)
        newFileGeno.close()
        newFilePheno.close()

else:
    while (phenoPosList==[] or phenoNegList==[]):
        alpha=random.uniform(-1,1)
        beta1=random.uniform(-1,1)
        beta2=random.uniform(-1,1)
        beta3=random.uniform(-1,1)
        beta12=random.uniform(-1,1)
        beta23=random.uniform(-1,1)
        beta13=random.uniform(-1,1)
        beta123=random.uniform(-1,1)

        for X1 in range (0,3):
            for X2 in range (0,3):
                for X3 in range (0,3):
                    value=alpha+beta1*X1+beta2*X2+beta3+X3+beta12*X1*X2+beta23*X2*X3+beta13*X1*X3+beta123*X1*X2*X3
                    phiList.append(value)
                    theMiddle=[X1,X2,X3,value]
                    globalList.append(theMiddle)



        for i in range(0,len(phiList)):
        	globalList[i][3]=1/(1+(math.exp(-phiList[i])))
        	print(globalList[i][3])
        	if (globalList[i][3]>0.5) :
        		phenoPosList.append([globalList[i][0],globalList[i][1],globalList[i][2]])
        	else :
        		phenoNegList.append([globalList[i][0],globalList[i][1],globalList[i][2]])

    for i in range(1,nbFiles+1):
        newFileGeno=open("{0}/{1}_Genotype_{2}.csv".format(outFolder,prefix,i),"w")
        header=""
        for j in range(1,nbVar-2):
            header=header+"N"+str(j)+","
        header=header+"M0P1,M0P2,M0P3\n"
        newFileGeno.write(header)
        countPosNeg=0
        incr=0
        for k in range(nbTotal):
            line=""
            for j in range(nbVar-3):
                line=line+str(random.randrange(3))+","
            if (countPosNeg<nbPositive):
                scramble=random.randrange(100)
                if (scramble<(scramblePercent)):
                    whichTrio=random.randrange(len(phenoNegList))
                    line=line+str(phenoNegList[whichTrio][0])+","+str(phenoNegList[whichTrio][1])+","+str(phenoNegList[whichTrio][2])
                    incr+=1
                else :
                    whichTrio=random.randrange(len(phenoPosList))
                    line=line+str(phenoPosList[whichTrio][0])+","+str(phenoPosList[whichTrio][1])+","+str(phenoPosList[whichTrio][2])
            else :
                if (incr>0):
                    whichTrio=random.randrange(len(phenoPosList))
                    line=line+str(phenoPosList[whichTrio][0])+","+str(phenoPosList[whichTrio][1])+","+str(phenoPosList[whichTrio][2])
                    incr-=1
                else :
                    whichTrio=random.randrange(len(phenoNegList))
                    line=line+str(phenoNegList[whichTrio][0])+","+str(phenoNegList[whichTrio][1])+","+str(phenoNegList[whichTrio][2])
            line=line+"\n"
            newFileGeno.write(line)
            countPosNeg+=1
        newFileGeno.close()

        newFilePheno=open("{0}/{1}_Phenotype_{2}.csv".format(outFolder,prefix,i),"w")
        newFilePheno.write("Class\n")
        for o in range(nbPositive):
            newFilePheno.write("1\n")
        for o in range(nbNegative):
            newFilePheno.write("0\n")
        newFilePheno.close()

print("Cas",phenoPosList)
print("Témoins",phenoNegList)
print(len(phenoNegList)+len(phenoPosList))

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
for i in range(1,nbFiles+1):
	newFile=open("{0}/{1}_Genome_{2}.csv".format(outFolder,prefix,i),"a")
	header=""
	for j in range(1,nbVar+1):
		header=header+"N"+str(j)+","
	header=header[0:-1]+"\n"
	line=""
	for k in range(1,nbTotal+1):
		for j in range(1,nbVar+1):
			line=line+str(random.randrange(3))+","
		line=line[0:-1]+'\n'
	newFile.write(header)
	newFile.write(line)

	newFile.close
alpha=random.uniform(-1,1)
beta1=random.uniform(-1,1)
beta2=random.uniform(-1,1)
beta12=random.uniform(-1,1)

phiList=[]
globalList=[]
for X1 in range (0,3):
	for X2 in range (0,3):
		value=alpha+beta1*X1+beta2*X2+beta12*X1*X2
		phiList.append(value)
		theMiddle=[X1,X2,value]
		globalList.append(theMiddle)


minim=min(phiList)
maxim=max(phiList)

phenoPosList=[]
phenoNegList=[]
for i in range(0,len(phiList)):
	globalList[i][2]=1/(1+(math.exp(-phiList[i])))
	print(globalList[i][2])
	if (globalList[i][2]>0.5) :
		phenoPosList.append([globalList[i][0],globalList[i][1]])
	else :
		phenoNegList.append([globalList[i][0],globalList[i][1]])

print("\n",alpha)
print(beta1)
print(beta2)
print(beta12)
print("Cas",phenoPosList)
print("Témoins",phenoNegList)

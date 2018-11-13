import re
import os
import random
import textwrap
import sys
import warnings

outFolder=input("Nom du répertoire de sortie :\n")
nbFiles=int(input("Nombre de fichier générés :\n"))
prefix=input("Préfixer commun aux fichiers générés :\n")
nbVar=int(input("Nombre des SNP à simuler :\n"))
nbPositive=int(input("Nombre de cas à simuler :\n"))
nbNegative=int(input("Nombre de témoins à simuler :\n"))
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

S1=random.uniform(-1,1)
S2=random.uniform(-1,1)
S12=random.uniform(-1,1)

phiList=[]
globalList=[]
for X1 in range (0,3):
	for X2 in range (0,3):
		value=S1*X1+S2*X2+S12*X1*X2
		phiList.append(value)
		theMiddle=[X1,X2,value]
		globalList.append(theMiddle)


minim=min(phiList)
maxim=max(phiList)

phenoPosList=[]
phenoNegList=[]
for i in range(0,len(phiList)):
	globalList[i][2]=(phiList[i]-minim)/(maxim-minim)

	if (globalList[i][2]>0.5) :
		phenoPosList.append([globalList[i][0],globalList[i][1]])
	else : 
		phenoNegList.append([globalList[i][0],globalList[i][1]])

# 1/(1+(-phi)²)

print(S1)
print(S2)
print(S12)
print("Cas",phenoPosList)
print("Témoins",phenoNegList)
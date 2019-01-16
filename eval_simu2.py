#!/usr/bin/env python3
# coding: utf-8

import re
import os
import random
import math
import sys
import warnings

entFolder=sys.argv[1]
sorFolder=sys.argv[2]
patternSize=int(sys.argv[3])
n_runs=sys.argv[4]

TP=0
FN=0
FP=0
countTP=0
countFP=0

for filename in os.listdir(entFolder):
    with open(os.path.join(entFolder, filename), 'r') as results:
        ligne=results.readline()
        while ligne:
            print(re.findall("M0P",ligne))
            if len(re.findall("M0P",ligne))==patternSize:
                countTP+=1
            else:
                countFP+=1
            ligne=results.readline()

eval=open(os.path.join(sorFolder, "eval_"+filename),'a')
if countTP!=0:
    TP+=1
    eval.write("TP"+'\n')
elif countFP!=0:
    countFP+=1
    eval.write("FP"+'\n')
else:
    FN+=1
    eval.write("FN"+'\n')
eval.close()


#recall = TP/(TP+FN)
#precision = TP/(TP+FP)
#fmeasure=2/(1/recall+1/precision)
#power = TP/n_runs
#power=open(os.path.join(sorFolder, "power_"+filename),'w')
#power.write("Recall="+str(recall)+'\n'+"Precision="+str(precision)+'\n'+"F-measure="+str(fmeasure)+'\n'+"Power="+str(power))
#power.close()

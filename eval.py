#!/usr/bin/env python3
# coding: utf-8

import re
import os
import random
import math
import sys
import warnings

entFolder=sys.argv[1]
entFile=sys.argv[2]
sorFolder=sys.argv[3]
patternSize=int(sys.argv[4])
causID=sys.argv[5]

TP=0
FN=0
FP=0
countTP=0
countFP=0


with open(os.path.join(entFolder,entFile), 'r') as results:
    ligne=results.readline()
    while ligne:
        if len(re.findall(str(causID),ligne))>=patternSize:
            countTP+=1
        else:
            countFP+=1
        if len(re.findall("seconds",ligne))>=0:
            timer=ligne
        ligne=results.readline()

eval=open(os.path.join(sorFolder, "results_"+entFile),'a')
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

time=open(os.path.join(sorFolder, "time"),'a')
time.write(timer)
time.close()

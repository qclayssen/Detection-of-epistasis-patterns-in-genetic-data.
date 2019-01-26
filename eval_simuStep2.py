#!/usr/bin/env python3
# coding: utf-8

import re
import os
import random
import math
import sys
import warnings

Folder=sys.argv[1]

TP=0
FN=0
FP=0
nb_runs=1
for filename in os.listdir(Folder):
    with open(os.path.join(Folder, filename), 'r') as results:
        ligne=results.readline()
        while ligne:
            if ligne.strip()=="TP":
                TP+=1
            elif ligne.strip()=="FP":
                FP+=1
            else:
                FN+=1
            nb_runs+=1
            ligne=results.readline()

if TP==0:
    recall=0
    precision=0
    fmeasure=0
    power=0
else:
    recall = TP/(TP+FN)
    precision = TP/(TP+FP)
    fmeasure=2/(1/recall+1/precision)
    power = TP/nb_runs

powerRes=open(os.path.join(Folder, "power_"+filename),'w')
powerRes.write("Recall="+str(recall)+'\n'+"Precision="+str(precision)+'\n'+"F-measure="+str(fmeasure)+'\n'+"Power="+str(power))
powerRes.close()

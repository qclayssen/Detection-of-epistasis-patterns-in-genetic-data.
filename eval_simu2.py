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
n_runs=sys.argv[3]
patternSize=sys.arg[4]

TP=0
FN=0
FP=0

for filename in os.listdir(entFolder):
    with open(os.path.join(entFolder, filename), 'r') as result:
        with open(os.path.join(repertoisortie, "eval_"+filename),'w') as eval:
            for ligne in result:
                if not ligne.strip():
                    FN+=1
                if re.findall("M0P",ligne):
                    TP+=1
                else:
                    FP+=1
            if FN!=0:





recall = TP/(TP+FN)
precision = TP/(TP+FP)
fmeasure =2/(1/recall+1/precision)
#power = TP/n_runs
power="1"
with open(os.path.join(repertoisortie, "power_"+filename),'w') as power:
    power.write("recall="+str(recall)+'\n'+"precision="+str(precision)+'\n'+"f-measure="+str(fmeasure)+'\n'+"power="+str(power))

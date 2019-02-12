#!/usr/bin/env python3
# coding: utf-8

import os
os.system("python3 simu_naive.py simupath3 1 simu2 50 1000 1000 3 5")

# 7 Arguments :
# python3 simu_naive.py {1} {2} {3} {4} {5} {6} {7}
# {1} = string : Output Folder
# {2} = int : Number of Files to generate
# {3} = string : Prefix of files
# {4} = int : Total Number of SNP
# {5} = int : Number of Cases (Phenotype 1)
# {6} = int : Number of Controls (Phenotype 0)
# {7} = int : Size of causal pattern (2 or 3)
# {8} = int : Percentage of scrambling (chance for a Control to get causal pattern)

#!/bin/bash

##############################

#This scripts takes 3 arguments :
#     {1} : Path to the Folder where are placed all data Files
#     {2} : Size of the pattern you're looking for (2 or 3)
#     {3} : The identifier of Causal Patterns in the genotype file (for the evaluation) ex:CAUS,M0P,X,...

##############################

#Example of command line :
# ./launch_meme.sh ../simu_Damien/Simu_naive_2snp_0.25 2 X


rm memetic
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results

make
simu=$(ls $1 | grep -i genotype)




for genotype in ${simu};
do
  for i in `seq 1 100`;
  do
    phenotype=$(echo ${genotype} | sed 's/Genotype/Phenotype/' | sed 's/genotype/phenotype/')
    ./memetic $1/${genotype} $1/${phenotype} #Execuction of the Method
    evalFile=$(ls ./outputs |grep ${genotype}) #First part of the evalutation : Creation of the 'results' file filed with TP/FP/FN
    ../eval.py outputs ${evalFile} results $2 $3
  done
  powerFile=$(ls ./results |grep ${genotype})
  ../eval_Step2.py results ${powerFile} #Second part of the evaluation : For each file
done

#!/bin/bash
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results


simu=$(ls data | grep -i genotype)

for genotype in ${simu};
do
  for i in `seq 1 10`;
    do
    phenotype=$(echo ${genotype} | sed 's/Genotype/Phenotype/' | sed 's/genotype/phenotype/')
    ../path_relinking data/${genotype} data/${phenotype} #Execuction of the Method
    evalFile=$(ls ./outputs |grep ${genotype}) #First part of the evalutation : Creation of the 'results' file filed with TP/FP/FN
    ../../eval.py outputs ${evalFile} results 2 M0P
    done
  powerFile=$(ls ./results |grep ${genotype})
  ../../eval_Step2.py results ${powerFile} #Second part of the evaluation : For each file
done

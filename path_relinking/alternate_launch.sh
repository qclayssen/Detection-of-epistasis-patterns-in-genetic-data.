#!/bin/bash
#rm path_relinking
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results

#make
simu=$(ls $1 | grep -i genotype)
genotypeFolder=$1
sizePattern=$2
identifierGen=$3

foo(){
  local run=$1
  local runtime=0
  for i in `seq 1 10`;
  do
    phenotype=$(echo $1 | sed 's/Genotype/Phenotype/' | sed 's/genotype/phenotype/')
    start=`date +%s`
    ./path_relinking ${genotypeFolder}/$1 ${genotypeFolder}/${phenotype} #Execuction of the Method
    end=`date +%s`
    diff=$((end-start))
    runtime=$(runtime+diff)
    evalFile=$(ls ./outputs |grep $1) #First part of the evalutation : Creation of the 'results' file filed with TP/FP/FN
    ../eval.py outputs ${evalFile} results $2 $3
  done
  powerFile=$(ls ./results |grep ${genotype})
  ../eval_Step2.py results ${powerFile} #Second part of the evaluation : For each file
}

for genotype in ${simu};
do
  foo "$genotype" "$2" "$3" &
done
wait

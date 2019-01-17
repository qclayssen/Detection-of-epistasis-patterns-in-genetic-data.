#!/bin/bash

rm memetic
rm -rf outputs

mkdir outputs

# Compilation if needed
make

# SMMB command line :
# ./SMMB <path_to_genotypes> <path_to_phenotypes>
#./memetic ./simupath3/simu2_Genotype_1.csv simupath3/simu2_Phenotype_1.csv
#../simu_Damien/Simu_naive_2snp_0.5/
simu=$(ls $1/*Genotype*)

for genotype in ${simu};
do
  for i in `seq 1 2`;
  do

  phenotype=$(echo ${genotype} | sed 's/Genotype/Phenotype/')
  ./memetic ${genotype} ${phenotype}

  python ../eval_simu.py outputs eval_${file}
  done
done

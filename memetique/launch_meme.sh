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
simu=$(ls $1 | grep -i genotype)

for genotype in ${simu};
do
  for i in `seq 1 10`;
  do

  phenotype=$(echo ${genotype} | sed 's/Genotype/Phenotype/')
  ./memetic ${genotype} ${phenotype}
  ../eval_simu.py outputs results 2
  #HEAPPROFILE=/tmp/netheap ./memetic ${genotype} ${phenotype}
  #HEAPCHECK=normal ./memetic ${genotype} ${phenotype}
  #valgrind --tool=callgrind --trace-children=yes ./memetic ${genotype} ${phenotype}

  done
done

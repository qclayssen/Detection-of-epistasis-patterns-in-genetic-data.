#!/bin/bash

rm memetic
rm -rf outputs

mkdir outputs

# Compilation if needed
make

# SMMB command line :
# ./SMMB <path_to_genotypes> <path_to_phenotypes>
#./memetic ./simupath3/simu2_Genotype_1.csv simupath3/simu2_Phenotype_1.csv
./memetic ./toy_dataset/genotypes_toy_dataset.txt ./toy_dataset/phenotypes_toy_dataset.txt

python ../eval_simu.py outputs eval_genotypes_toy_dataset

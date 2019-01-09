#!/bin/bash

rm memetic

# Compilation if needed
make

# SMMB command line :
# ./SMMB <path_to_genotypes> <path_to_phenotypes>
#./memetic ./simupath3/simu2_Genotype_1.csv simupath3/simu2_Phenotype_1.csv
./memetic ./toy_dataset/genotypes_toy_dataset.txt ./toy_dataset/phenotypes_toy_dataset.txt

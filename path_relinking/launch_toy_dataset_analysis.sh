#!/bin/bash
rm path_relinking
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results
# Compilation if needed
make

# Path_relinking command line :
# ./àath_relinking <path_to_genotypes> <path_to_phenotypes>
./path_relinking ./toy_dataset/simu3_Genotype_1.csv ./toy_dataset/simu3_Phenotype_1.csv

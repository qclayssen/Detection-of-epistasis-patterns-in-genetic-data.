#!/bin/bash
rm memetic
rm -rf outputs

mkdir outputs
# Compilation if needed
make

# Path_relinking command line :
# ./àath_relinking <path_to_genotypes> <path_to_phenotypes>
./memetic ./genotype_model1_1_01p_0005h_01m_002.txt ./phenotype_model1_1_01p_0005h_01m_002.txt
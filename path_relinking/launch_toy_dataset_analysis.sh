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
./path_relinking ./toy_dataset/genotypes_toy_dataset.txt ./toy_dataset/phenotypes_toy_dataset.txt

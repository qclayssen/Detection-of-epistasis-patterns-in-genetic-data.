#!/bin/bash
rm memetic
rm -rf outputs

mkdir outputs
# Compilation if needed
make

# Path_relinking command line :
# ./àath_relinking <path_to_genotypes> <path_to_phenotypes>
./memetic ./toy_dataset/genotypes_toy_dataset.txt ./toy_dataset/phenotypes_toy_dataset.txt

#!/bin/bash

# Compilation if needed
make

# SMMB command line :
# ./SMMB <path_to_genotypes> <path_to_phenotypes>
./path_relinking ../simu2/simu2_Genotype_1.csv ../simu2/simu2_Phenotype_1.csv

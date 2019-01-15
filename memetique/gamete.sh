#!/bin/bash


b=$(basename $1)
mkdir file_$1
awk '{print $NF}' $1 > ./file_$1/phenotype_${b}
awk '{$NF=""; print $0}' $1 | sed 's/[[:space:]]/,/g'> ./file_$1/genotype_${b}

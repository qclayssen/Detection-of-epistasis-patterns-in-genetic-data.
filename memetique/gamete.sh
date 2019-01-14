#!/bin/bash


b=$(basename $1)
awk '{print $NF}' $1 > phenotype_${b}
awk '{$NF=""; print $0}' $1 > genotype_${b}

#!/bin/bash

#./gamete.sh ./repository/model1/model1_0_01p_0005h_005m/

b=$(basename $1)
mkdir data_gametes_${b}
model=$(ls $1)
for file in $model
do
awk '{print $NF}' $1/${file} > ./data_gametes/phenotype_${file}
awk '{$NF=""; print $0}' $1/${file} | sed 's/[[:space:]]/,/g' | sed 's/,$//' > ./data_gametes_${b}/genotype_${file}
done

rm path_relinking
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results

make
for i in `seq 1 100`;
do
  ./path_relinking ./toy_dataset/simu3_Genotype_1.csv ./toy_dataset/simu3_Phenotype_1.csv
  ../eval_simu2.py outputs results 2 100
done

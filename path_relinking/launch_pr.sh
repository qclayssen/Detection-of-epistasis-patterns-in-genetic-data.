rm path_relinking
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results

make
for i in `seq 1 10`;
do
  ./path_relinking ./toy_dataset/simu3_Genotype_1.csv ./toy_dataset/simu3_Phenotype_1.csv
  ../eval_simu.py outputs results 3
done
../eval_simuStep2.py results

rm path_relinking
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results

make
for i in `seq 1 10`;
do
  ./path_relinking ./toy_dataset/genotype_model1_1_01p_0005h_01m_002.txt ./toy_dataset/phenotype_model1_1_01p_0005h_01m_002.txt
  ../eval_simu.py outputs results 2
done
../eval_simuStep2.py results

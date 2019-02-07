#rm path_relinking
rm -rf outputs
rm -rf results

mkdir outputs
mkdir results

#make
simu=$(ls $1 | grep -i genotype)

for genotype in ${simu};
do
  for i in `seq 1 10`;
  do
    phenotype=$(echo ${genotype} | sed 's/Genotype/Phenotype/' | sed 's/genotype/phenotype/')
    ./path_relinking ./toy_dataset/${genotype} ./toy_dataset/${phenotype}
    ../eval_simu.py outputs results 2
  done
  ../eval_simuStep2.py results
done

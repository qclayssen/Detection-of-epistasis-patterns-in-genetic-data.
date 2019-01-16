make
for i in `seq 1 100`;
do
  ./path_relinking ./toy_dataset/genotypes_toy_dataset.txt ./toy_dataset/phenotypes_toy_dataset.txt
done

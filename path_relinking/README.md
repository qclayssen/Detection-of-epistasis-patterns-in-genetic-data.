This path-relinking algorithm was realized during a Master 2 project, and it's purpose was to be applied to epistasy, and to judge of it's efficiency on this kind of datas.

*Main contributors :*
Quentin Clayssen and Antoine Lainé, Master 2 students at University of Nantes, in Bioinformatics.

*Other contributors :*
Christine Sinoquet (professor, and supervisor of the project)
Clément Niel (part of SMMB-ACO script where used as a basis)



**Tutorial :**

*1-Installation*
No step of installation required. The only thing you need to know is the path where the Folder was downloaded/placed.
Every C++/Python library used to run the method should be installed natively.

Note : The Boost library is directly in the archive, so it doesn't need to be installed.

*2-Compilation*
Be sure you are in the path-relinking directory. Here you should be able to use the command :
    make

*3-Execution*
The method can be launched simply by using the command :
    ./launch_pr.sh genotype_file.txt phenotype_file.txt
Replacing the parameters with the path to your genotype and phenotype files.

Exectution on toy example :
    ./launch_pr.sh toy_example/data/genotype_model1_0_01p_0005h_005m_001.txt toy_example/data/phenotype_model1_0_01p_0005h_005m_001.txt

*4-Results*
The results can be found in the 'outputs' folder, under this format :
==========================================
Pattern	p-value	score
<N24,N26,N32>	8.67208e-05	62.1096
<N23,N35,N93>	0.000739019	55.0874
<N6,N35,N93>	0.00230297	51.1236
<N5,N15,N23>	0.00824289	46.3923
<N17,N26,N73>	0.0234011	42.2015

durée: 0.834144 seconds
==========================================

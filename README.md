# Detection-of-epistasis-patterns-in-genetic-data.
Master2 project Bioinformatics, University of Nantes
By Quentin Clayssen, Antoine Laine 

Detection of epistasis patterns using [Memetic algorithm](https://en.wikipedia.org/wiki/Memetic_algorithm) and Path Relinking method



## Initialisation of project

In order to compile the whole project a makefile is provided at project's root. This makefile will call each method's makefile and produce both executables. This makefile, can recreate directory structure if needed.


### Compilation instruction
To recreate directory structure for both methods, please call this line at the root of the project:

    make install

To compile both methods, please call this lines at the root of the project:

    make

If a recompilation is needed please use to purge all compiled files:

    make clean


## Generate a naive data set

This project is provided with a tool to generate some naive data set. This tool fit with a logistic regression and produce a data set of given sizes with a given number of causal SNPs. This operation is performed by simu_naive.py.



    ./simu_naive.py -p simu_naive -f 1 -v 28 -pa 2000 -o toy_dataset -c 2000 -s 2

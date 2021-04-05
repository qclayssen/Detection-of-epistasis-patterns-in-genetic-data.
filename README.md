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

### Prerequisite, generation of right virtual environment

These steps are done if installation_of_project.sh was used.

Some packages are required to run this tool. List of them is provided with needed version in environment.yml or environment.txt file. In order to generate a virtual environment where simu_naive.py can be executed please follow next steps.

For this step conda is needed, please follow instructions from anaconda documentation to get it installed. If conda is installed please execute this line to generate an environment called projet_c:

    conda env create -f ./environment_to_execute_python/environment.yml

If virtualenv if prefered and is not installed. Please use following lines:

    pip install virtualenv
    virtualenv -p /usr/bin/python3.6 projet_c
    source projet_c/bin/activate
    pip install -r environment.txt
    deactivate



Following this, activate this environment and execute simu_naive with arguments:

    source activate projet_c
    ./simu_naive.py -p simu_naive -f 1 -v 28 -pa 2000 -o toy_dataset -c 2000 -s 2
    source deactivate

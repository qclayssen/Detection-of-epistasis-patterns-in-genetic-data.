BOOST_FOLDER=/home/antoine/Documents/project_c_sinoquet/include/boost_1_68_0

SRC=path_relinking.cpp\
../src/stats.cpp\
../src/tools.cpp\
../src/path_relinking_func.cpp

FLAG=-Wall -std=c++11 -I$(BOOST_FOLDER)

all: ./path_relinking

path_relinking : $(SRC)
	@g++ $^ -o $@  $(FLAG)

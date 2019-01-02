BOOST_FOLDER=../include/boost_1_68_0


SRC=memetic_funk.cpp\
../src/Parameters_file_parsing.cpp\
../src/stats.cpp\
../src/tools.cpp\
../src/path_relinking_func.cpp\
../src/memetic_funk.cpp

FLAG=-Wall -std=c++11 -I$(BOOST_FOLDER)

all: ./path_relinking

path_relinking : $(SRC)
	@g++ $^ -o $@  $(FLAG)

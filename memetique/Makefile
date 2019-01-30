BOOST_FOLDER=../include/boost_1_68_0


SRC=memetic.cpp\
../src/Parameters_file_parsing.cpp\
../src/stats.cpp\
../src/tools.cpp\
../src/path_relinking_func.cpp\
../src/memetic_funk.cpp

FLAG=-Wall -std=c++11 -O2 -I$(BOOST_FOLDER)

all: ./memetic

memetic : $(SRC)
	@g++ $^ -o $@  $(FLAG)

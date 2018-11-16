all: parsingmake

INC=-I$(BOOST_FOLDER) \
 -I./include \
 -I./include/parsing \
 -I./include/services \
 -I./include/structures \
 -I./include/statistics \
 -I./include/information_theory \
 -I./include/usecases \
 -I./include/markov


parsingmake: parameters_file.o
     g++ -c parameters_file_parsing.cpp c++11 -fopenmp -Wall -Wextra -DNDEBUG -g $(INC)

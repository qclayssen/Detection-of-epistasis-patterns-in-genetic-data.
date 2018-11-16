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
	g++ -o parameters_file_parsing.o -c parameters_file_parsing.cpp -std=c++11 -Wall -Wextra -DNDEBUG
target: parameters_file.o
     g++ -o exec_parser parameters_file_parsing.o -std=c++11 -Wall -Wextra -DNDEBUG

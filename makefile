all: exec_parser

exec_parser: parameters_file.o main.o
	g++ -o exec_parser parameters_file_parsing.o main.o -std=c++11 -Wall -Wextra -DNDEBUG

parameters_file.o: parameters_file_parsing.cpp
	g++ -o parameters_file_parsing.o -c parameters_file_parsing.cpp -std=c++11 -Wall -Wextra -DNDEBUG

main.o: main.cpp
	g++ -o main.o -c main.cpp -std=c++11 -Wall -Wextra -DNDEBUG

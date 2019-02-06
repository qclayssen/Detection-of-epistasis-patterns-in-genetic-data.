#ifndef PATH_RELINKING_FUNC_HPP
#define PATH_RELINKING_FUNC_HPP

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <list>
#include <numeric>
#include <random>
#include <vector>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/algorithm/string.hpp>


#include "typesandstruct.hpp"
#include "stats.hpp"
#include "tools.hpp"
#include "CSVParser.hpp"
#include "common.h"

using namespace std;
using namespace std::chrono;

vector<patternscore> initialize_elite_solutions(unsigned int k,vector<patternscore> patternscoreList);

vector<patternscore> select_two_solutions_at_random(vector<patternscore> elite_sols);

patternscore select_closest_neighbor_to_guiding_solution(patternscore s,patternscore sB, vector<patternscore> patternscoreList, int s_n);

int promizing_score(patternscore s_closest_neighbour,vector<patternscore> elite_sols);

void update(patternscore s_opt, vector<patternscore>* adr_elite_sols);



#endif

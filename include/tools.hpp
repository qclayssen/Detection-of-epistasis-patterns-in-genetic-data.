#ifndef TOOLS_HPP
#define TOOLS_HPP

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
#include <libgen.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/algorithm/string.hpp>

#include "typesandstruct.hpp"
#include "stats.hpp"
#include "CSVParser.hpp"
#include "common.h"

using namespace std;
using namespace std::chrono;

vector<string> get_snp_list(string genos_file);

void cout_list(vector<patternscore> list_to_cout,vector<string> snpNameList);

int calculate_delta(patternscore s, patternscore sB);

vector<patternscore> neighbours(patternscore s,vector<patternscore> patternscoreList);

patternscore hill_climbing_lc(patternscore s_closest_neighbour, vector<patternscore> patternscoreList,blas_matrix genos,blas_matrix phenos_m);

void outfile(string genos_file,vector<string> snpNameList,vector<patternscore> best_solutions);


#endif

#ifndef MEMETIC_FUNK_HPP
#define MEMETIC_FUNK_HPP

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

vector<parents_pairs> select_pairs_of_individuals_to_be_crossed(vector<patternscore> n_pairs_selected_parents);


void perform_one_mutation_per_child(vector<patternscore>* adr_children_parents,int prob_mutation);

vector<patternscore> create_two_children_for_each_selected_pair_of_parents(vector<parents_pairs> pairs_of_parents);

float add_gtest_score (patternscore pattern,blas_matrix genos,blas_matrix phenos_m);

void update(patternscore s_opt, vector<patternscore>* adr_elite_sols);

vector<patternscore> identify_best_solutions(vector<patternscore> pop, int k, int n);

vector<patternscore> initialize_population(int n,vector<patternscore> patternscoreList);

void update_population(vector<patternscore> children_parents, vector<patternscore>* adr_pop,int n);

patternscore hill_climbing_lc2(patternscore s_closest_neighbour, vector<patternscore> patternscoreList);

bool compareByLength(const patternscore &a, const patternscore &b);


#endif

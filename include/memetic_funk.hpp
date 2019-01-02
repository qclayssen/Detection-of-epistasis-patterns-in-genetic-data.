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

float add_gtest_score (patternscore pattern,blas_matrix genos,blas_matrix phenos_m)

void update(patternscore s_opt, vector<patternscore>* adr_elite_sols)

vector<patternscore> identify_best_solutions(vector<patternscore> pop, int k, int n)

vector<patternscore> initialize_population(int n,vector<patternscore> patternscoreList)

void update_population(vector<patternscore> children_parents, vector<patternscore>* adr_pop,int n)

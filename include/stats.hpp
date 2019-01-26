#ifndef STATS_HPP
#define STATS_HPP

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
#include "CSVParser.hpp"
#include "common.h"

using namespace std;
using namespace std::chrono;

void create_contingency_table_pattern2(int l1,int l2,contingence2SNP* adr_contingence,blas_matrix genos, blas_matrix phenos_m);

void create_contingency_table_pattern3(int l1,int l2,int l3,contingence3SNP* adr_contingence,blas_matrix genos, blas_matrix phenos_m);

score_pval g_test_2SNP(contingence2SNP contingence2);

score_pval g_test_3SNP(contingence3SNP contingence2);

float add_gtest_pval(patternscore pattern,blas_matrix genos,blas_matrix phenos_m);

float add_gtest_score(patternscore pattern,blas_matrix genos,blas_matrix phenos_m);

score_pval add_gtest_results (patternscore pattern,blas_matrix genos,blas_matrix phenos_m);

#endif

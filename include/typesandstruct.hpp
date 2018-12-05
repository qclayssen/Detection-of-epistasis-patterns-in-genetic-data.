#ifndef TYPES_AND_STRUCT_HPP
#define TYPES_AND_STRUCT_HPP

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <numeric>

#include "common.h"

using namespace std;
using namespace std::chrono;

typedef int contingence2SNP[3][10];
typedef int contingence3SNP[3][28];
struct patternscore {
  int snp1;
  int snp2;
  int snp3;
  string pattern1;
  string pattern2;
  string pattern3;
  float score = -1;
};

#endif

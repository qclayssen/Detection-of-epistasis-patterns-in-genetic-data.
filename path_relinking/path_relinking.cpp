//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
//Parsing of CSV Files for Epistasis detection
//Created :09/11/18
//Modified :11/02/2019

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

#include "../include/Parameters_file_parsing.hpp"
#include "../include/typesandstruct.hpp"
#include "../include/CSVParser.hpp"
#include "../include/stats.hpp"
#include "../include/tools.hpp"
#include "../include/path_relinking_func.hpp"
#include "../include/common.h"

using namespace std;
using namespace std::chrono;
typedef std::chrono::high_resolution_clock Clock;


int main(int argc, char *argv[])
{
  auto tinit = Clock::now(); //Initialize the timer
  if(argc < 3) //Verify that 2 files where given as parameters
  {
      cerr << "Missing parameter :\n"
           << "\t./path_relinking <path_to_genotypes> <path_to_phenotypes>"
           << endl;
      exit(-1);
  }

    // Arguments / Parameters set
    string genos_file = argv[1];
    string phenos_file = argv[2];

    parameters_file_parsing params;
    int header = params.header;
    char separator = params.separator;
    int s_n = params.s_n;
    unsigned int k = params.k;

    // Data importation, using boost (part of it was inspired by Clément Niel's version of SMMB-ACO)
    CSVParser<int> genos_csv(genos_file, separator, header);
    CSVParser<int> phenos_csv(phenos_file, separator, header);
    blas_matrix genos = genos_csv.data();
    blas_matrix phenos_m = phenos_csv.data();
    blas_column phenos(phenos_m, 0);

    // Creation of global parameters
    int l1,l2,l3;
    vector<string> snpNameList; //This vector contains all the names of the SNP (from the header of genotype file)
    snpNameList = get_snp_list(genos_file);
    vector<patternscore>patternscoreList;
    vector<patternscore>::iterator iterpatternscoreList;

    // Here is the main tool we use in our method : a vector of structure, with each structure corresponding to a solution(pattern).
    //The structure is names patternscore, and the vector is pattentscoreList
    //For now, we don't touch the score/p-value, as we only need to calculate them on solutions that we will use.
    for (l1=0;l1<int(genos.size2())-1;l1++){ //First SNP of the pattern
      for (l2=l1+1;l2<int(genos.size2());l2++){ //Second SNP of the pattern
        patternscore p1;
        p1.pattern1=snpNameList[l1];
        p1.pattern2=snpNameList[l2];
        p1.pattern3="";
        p1.snp1=l1;
        p1.snp2=l2;
        p1.snp3=-1;
        patternscoreList.push_back(p1);
      }
    }

    //In the same vector, we push the possible patterns of size 3.
    for (l1=0;l1<int(genos.size2())-2;l1++){ //First SNP of the pattern
      for (l2=l1+1;l2<int(genos.size2())-1;l2++){ //Second SNP of the pattern
        for(l3=l2+1;l3<int(genos.size2());l3++){ //Third SNP of pattern
          patternscore p1;
          p1.snp1=l1;
          p1.snp2=l2;
          p1.snp3=l3;
          p1.pattern1=snpNameList[l1];
          p1.pattern2=snpNameList[l2];
          p1.pattern3=snpNameList[l3];
          patternscoreList.push_back(p1);
        }
      }
    }


    score_pval biScore;
    vector<patternscore> elite_sols;
    elite_sols = initialize_elite_solutions(k,patternscoreList); //C.F path_relinking_func.cpp

    for (unsigned int i=0;i<elite_sols.size();i++){ //We set the score and p-value of each elite solution (through contingency and gtest). A small local search was added here to make the starting elite solutions a little more accurate.
        biScore=add_gtest_results(elite_sols[i],genos,phenos_m);
        elite_sols[i].score=biScore.score;
        elite_sols[i].pval=biScore.pval;
        elite_sols[i]=hill_climbing_lc(elite_sols[i],patternscoreList,genos,phenos_m,s_n);
    }
    vector<patternscore>* adr_elite_sols = &elite_sols; //Getting the adress were our elite solutions are stocked, so that we can modify them in a function.


    vector<patternscore> sA_sB;
    sA_sB=select_two_solutions_at_random(elite_sols); //C.F path_relinking_func.cpp

    patternscore s = sA_sB[0]; //sA is now named s
    patternscore sB = sA_sB[1]; //sB stays sB

    while (calculate_delta(s,sB)>0){
      //While s and sB are different, the neighbour of s that is the closest to sB is compared to all of our elite solutions.
      //If the neighbour is better than the lowest value, a local search is used.
      //The solution outputted of local search replaces the lowest elite solution.
      //Even if the elite solutions stay the same, s is replaced by the neighbour, until s equal sB.
      patternscore s_closest_neighbour=select_closest_neighbor_to_guiding_solution(s,sB,patternscoreList, patternscoreList.size()); //C.F path_relinking_func.cpp
      biScore=add_gtest_results(s_closest_neighbour,genos,phenos_m);
      s_closest_neighbour.score=biScore.score;
      s_closest_neighbour.pval=biScore.pval;
      if (promizing_score(s_closest_neighbour,elite_sols)==1){
        //cout<<"Local search"<<endl;
        patternscore s_opt=hill_climbing_lc(s_closest_neighbour,patternscoreList,genos,phenos_m,s_n);
        update(s_opt,adr_elite_sols);
      }
      else{
        //cout<<"No local search"<<endl;
      }
      s=s_closest_neighbour;

    }

    elite_sols=sort_solutions(elite_sols); //Sorts solution by p-value, or by score if p-value=0.

    auto tend = Clock::now(); //Checks timer, and calculates the difference between initial timer and actual timer
    float duree = duration_cast<duration<double>>(tend - tinit).count();
    cout << "Total execution time: "<<duree<< " seconds"<<endl;

    outfilePR(genos_file,snpNameList,elite_sols,s_n,duree); //Prints all the results in a file

    return 0;
}

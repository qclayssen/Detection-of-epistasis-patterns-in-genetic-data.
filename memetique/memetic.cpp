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
#include <stdio.h>
#include <ctime>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/algorithm/string.hpp>

#include "../include/Parameters_file_parsing.hpp"
#include "../include/typesandstruct.hpp"
#include "../include/CSVParser.hpp"
#include "../include/stats.hpp"
#include "../include/tools.hpp"
#include "../include/memetic_funk.hpp"
#include "../include/common.h"

using namespace std;
using namespace std::chrono;
typedef std::chrono::high_resolution_clock Clock;



int main(int argc, char *argv[])
{
  auto t1 = Clock::now();//Initialize the timer
  if(argc < 3) //Verify that 2 files where given as parameters
  {
      cerr << "Missing parameter :\n"
           << "\t./memetic <path_to_genotypes> <path_to_phenotypes>"
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
    int n_p = params.n_pairs_selected_parents;
    unsigned int best_k = params.best_k;
    int prob_mutation = params.prob_mutation;// probability of mutation
    int n_it = params.n_it; // Nombre itération
    int n = params.n; // Size Population initiale
    int h = 0;
    score_pval biScore; //Holds the result of g-tests




  //  DATA IMPORTATION
    CSVParser<int> genos_csv(genos_file, separator, header);
    CSVParser<int> phenos_csv(phenos_file, separator, header);
    blas_matrix genos = genos_csv.data();
    blas_matrix phenos_m = phenos_csv.data();
    blas_column phenos(phenos_m, 0);
    int l1,l2,l3;



    vector<string> snpNameList;
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




    vector<patternscore>* adr_patternscoreList = &patternscoreList;

    //create population
    vector<patternscore>pop=initialize_population(n, adr_patternscoreList);;

    // add score/pval to individues of population
    vector<patternscore>* adr_pop = &pop;
    for (int i=0;i<pop.size();i++){
        biScore=add_gtest_results(pop[i],genos,phenos_m);
        pop[i].score=biScore.score;
        pop[i].pval=biScore.pval;
    }


    //local search for each individues of population
    for (int l=0;l<pop.size();l++){
            patternscore s_opt=hill_climbing_lc(pop[l],patternscoreList,genos,phenos_m,s_n);
            pop[l]=s_opt;
          }

//stop condition for evolution of population
    while (h < n_it){

      vector<patternscore> pop_init=pop;
      for (int l=0;l<pop.size();l++){
              pop[l].idparent=l;
            }


      vector<parents_pairs> pairs_of_parents=select_pairs_of_individuals_to_be_crossed(pop,n_p);

      vector<patternscore> children_parents;
      children_parents=create_two_children_for_each_selected_pair_of_parents(pairs_of_parents);
      vector<patternscore>* adr_children_parents = &children_parents;
      perform_one_mutation_per_child(adr_children_parents,prob_mutation,patternscoreList);

      for (int i=0;i<children_parents.size();i++){
        if (children_parents[i].score==0 && children_parents[i].pval==0)
        {biScore=add_gtest_results(children_parents[i],genos,phenos_m);
        children_parents[i].score=biScore.score;
        children_parents[i].pval=biScore.pval;}
      }

      //replace parent by children with higher score
      update_population(children_parents, adr_pop,n);


      //local search on new population
      for (int o=0;o<pop.size();o++){
              patternscore s_opt=hill_climbing_lc(pop[o],patternscoreList,genos,phenos_m,s_n);
              pop[o]=s_opt;
            }
      h=h+1;
    }
    vector<patternscore> best_solutions = identify_best_solutions(pop,best_k,n);
    auto t5 = Clock::now();
    float duree =duration_cast<duration<double>>(t5 - t1).count();
    //Checks timer, and calculates the difference between initial timer and actual timer

    outfile(genos_file,snpNameList, best_solutions,s_n,n,duree,n_it);//Prints all the results in a file

    return 0;



    }

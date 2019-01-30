//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
//Parsing of CSV Files for Epistasis detection
//Created :
//Modified :

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
#include "../include/path_relinking_func.hpp"
#include "../include/memetic_funk.hpp"
#include "../include/common.h"

using namespace std;
using namespace std::chrono;



int main(int argc, char *argv[])
{
  if(argc < 3)
  {
      cerr << "Missing parameter :\n"
           << "\t./path_relinking <path_to_genotypes> <path_to_phenotypes>"
           << endl;
      exit(-1);
  }

    // Arguments
    string genos_file = argv[1];
    string phenos_file = argv[2];

    parameters_file_parsing params;
    int header = params.header;
    char separator = params.separator;


  //  DATA IMPORTATION
    CSVParser<int> genos_csv(genos_file, separator, header);
    CSVParser<int> phenos_csv(phenos_file, separator, header);
    blas_matrix genos = genos_csv.data();
    blas_matrix phenos_m = phenos_csv.data();
    blas_column phenos(phenos_m, 0);
    //cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;
    int l1,l2,l3;



    vector<string> snpNameList;
    snpNameList = get_snp_list(genos_file);

    vector<patternscore>patternscoreList;
    vector<patternscore>::iterator iterpatternscoreList;




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
    int n_c = params.n_c;
    unsigned int k = params.k;
    int prob_mutation = params.prob_mutation;
    int n_it = params.n_it; // Nombre itération
    int n = params.n; // Size Population initiale
    int h = 1;
    score_pval biScore; //Holds the result of g-tests

    vector<patternscore>pop=initialize_population(n, patternscoreList);;
    cout<<"pop initiale:"<<endl;
    //cout_list(pop,snpNameList);
    vector<patternscore>* adr_pop = &pop;
    for (int i=0;i<pop.size();i++){
        biScore=add_gtest_results(pop[i],genos,phenos_m);
        pop[i].score=biScore.score;
        pop[i].pval=biScore.pval;
    }
cout<<"1"<<endl;
    for (int l=0;l<pop.size();l++){
            patternscore s_opt=hill_climbing_lc(pop[l],patternscoreList,genos,phenos_m,n_c);
            //cout<<s_opt.snp1<<endl;
            update(s_opt,adr_pop);
          }

      cout<<"pop après recherche:"<<endl;
      //std::time_t result = std::time(nullptr);
      //std::cout << std::asctime(std::localtime(&result))<<endl;
      //cout_list(pop,snpNameList);

    vector<patternscore>n_pairs_selected_parents=pop;
    for (int i=0;i<n_pairs_selected_parents.size();i++){
      biScore=add_gtest_results(n_pairs_selected_parents[i],genos,phenos_m);
      n_pairs_selected_parents[i].score=biScore.score;
      n_pairs_selected_parents[i].pval=biScore.pval;
    }
    while (h < n_it){
      vector<parents_pairs> pairs_of_parents=select_pairs_of_individuals_to_be_crossed(n_pairs_selected_parents);
      //cout_list_indiv(pairs_of_parents);

      vector<patternscore> children_parents;
      children_parents=create_two_children_for_each_selected_pair_of_parents(pairs_of_parents);
      //cout_list(children_parents);


      cout<<"enfant:"<<endl;
      //std::time_t result1 = std::time(nullptr);
      //std::cout << std::asctime(std::localtime(&result1))<<endl;
      //cout_list(children_parents,snpNameList);
      vector<patternscore>* adr_children_parents = &children_parents;
      perform_one_mutation_per_child(adr_children_parents,prob_mutation);

      for (int i=0;i<children_parents.size();i++){
        biScore=add_gtest_results(children_parents[i],genos,phenos_m);
        children_parents[i].score=biScore.score;
        children_parents[i].pval=biScore.pval;
      }
      cout<<"enfant muté:"<<endl;
      //std::time_t result2 = std::time(nullptr);
      //std::cout << std::asctime(std::localtime(&result2))<<endl;
      //cout_list(children_parents,snpNameList);
      //update_population(children_parents, adr_pop,n);
    //  cout<<"pop size:"<<pop.size()<<endl;
      for (int o=0;o<pop.size();o++){
              patternscore s_opt=hill_climbing_lc(pop[o],patternscoreList,genos,phenos_m,n_c);
              update(s_opt,adr_pop);
            }

      //cout<<"pop finale:"<<endl;
      //cout_list(pop,snpNameList);

      h=h+1;
    }
    //char filename= "out.txt";
    vector<patternscore> best_solutions = identify_best_solutions(pop,k,n);
    cout<<"pop finale trié:"<<endl;
    //cout_list(best_solutions,snpNameList);
    //std::time_t result3 = std::time(nullptr);
    //std::cout << std::asctime(std::localtime(&result3))<<endl;
    outfile(genos_file,snpNameList, best_solutions);

    return 0;



    }

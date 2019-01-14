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
    cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;
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
            p1.snp3=NULL;
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
            //  patternscoreList.push_back(p1);
            }
          }
        }


    unsigned int k = params.k;
    int prob_mutation = params.prob_mutation;
    int n_it = params.n_it; // Nombre itération
    int n = params.n; // Size Population initiale
    int h = 1;


    vector<patternscore>pop=initialize_population(n, patternscoreList);;
    cout<<"pop initiale:"<<endl;
    cout_list(pop,snpNameList);
    vector<patternscore>* adr_pop = &pop;
    for (int i=0;i<pop.size();i++){
        pop[i].score=add_gtest_score(pop[i],genos,phenos_m);
    }

    for (int l=0;l<pop.size();l++){
            patternscore s_opt=hill_climbing_lc(pop[l],patternscoreList,genos,phenos_m);
            //cout<<s_opt.snp1<<endl;
            update(s_opt,adr_pop);
          }

      cout<<"pop après recherche:"<<endl;
      cout_list(pop,snpNameList);
    vector<patternscore>n_pairs_selected_parents=pop;
    for (int i=0;i<n_pairs_selected_parents.size();i++){
        n_pairs_selected_parents[i].score=add_gtest_score(n_pairs_selected_parents[i],genos,phenos_m);
    }
    while (h < n_it){
      vector<parents_pairs> pairs_of_parents=select_pairs_of_individuals_to_be_crossed(n_pairs_selected_parents);
      //cout_list_indiv(pairs_of_parents);

      vector<patternscore> children_parents;
      children_parents=create_two_children_for_each_selected_pair_of_parents(pairs_of_parents);
      //cout_list(children_parents);


      for (int j=0;j<children_parents.size();j++){
        //  cout<<children_parents.size()<<"i"<<endl;
          float score=add_gtest_score(children_parents[j],genos,phenos_m);
      }
      cout<<"enfant:"<<endl;
      cout_list(children_parents,snpNameList);
      vector<patternscore>* adr_children_parents = &children_parents;
      perform_one_mutation_per_child(adr_children_parents,prob_mutation);

      for (int i=0;i<children_parents.size();i++){
          children_parents[i].score=add_gtest_score(children_parents[i],genos,phenos_m);
      }
      cout<<"enfant muté:"<<endl;
      cout_list(children_parents,snpNameList);
      update_population(children_parents, adr_pop,n);
      cout<<"pop size:"<<pop.size()<<endl;
      for (int o=0;o<pop.size();o++){
              patternscore s_opt=hill_climbing_lc(pop[o],patternscoreList,genos,phenos_m);
              update(s_opt,adr_pop);
            }

      cout<<"pop finale:"<<endl;
      cout_list(pop,snpNameList);

      h=h+1;
    }
    //char filename= "out.txt";
    vector<patternscore> best_solutions = identify_best_solutions(pop,k,n);
    cout<<"pop finale trié:"<<endl;
    cout_list(best_solutions,snpNameList);
    outfile( snpNameList, best_solutions);
    /*ofstream file( filename.c_str() );
    file << cout_list(best_solutions,snpNameList);
    file.close();


      if ((freopen(filename, "w", stdout)) != NULL)
  {
      cout_list(best_solutions,snpNameList);

      fclose (stdout);
  }
  else
  {
      cout<<"fail"<<endl;
  }*/
    return 0;



    }

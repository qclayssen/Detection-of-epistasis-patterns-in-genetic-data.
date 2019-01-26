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


    unsigned int k = params.k;
    vector<patternscore> elite_sols;
    elite_sols = initialize_elite_solutions(k,patternscoreList);
    for (unsigned int i=0;i<elite_sols.size();i++){
        elite_sols[i].pval=add_gtest_pval(elite_sols[i],genos,phenos_m);
        elite_sols[i]=hill_climbing_lc(elite_sols[i],patternscoreList,genos,phenos_m);
    }
    vector<patternscore>* adr_elite_sols = &elite_sols;
    cout<<"Elite solutions :"<<endl;
    cout_list(elite_sols,snpNameList);

    vector<patternscore> sA_sB;
    sA_sB=select_two_solutions_at_random(elite_sols);
    cout<<"Two random : (sA, sB)"<<endl;
    cout_list(sA_sB,snpNameList);
    patternscore s = sA_sB[0];
    patternscore sB = sA_sB[1];

    while (calculate_delta(s,sB)>0){
      patternscore s_closest_neighbour=select_closest_neighbor_to_guiding_solution(s,sB,patternscoreList);
      s_closest_neighbour.pval=add_gtest_pval(s_closest_neighbour,genos,phenos_m);
      if (promizing_score(s_closest_neighbour,elite_sols)==1){
        cout<<"Recherche locale"<<endl;
        patternscore s_opt=hill_climbing_lc(s_closest_neighbour,patternscoreList,genos,phenos_m);
        update(s_opt,adr_elite_sols);
      }
      else{
        cout<<"Pas de recherche locale"<<endl;
      }
      s=s_closest_neighbour;
    }
    cout<<endl<<"Solutions d'élite finales:"<<endl;
    cout_list(elite_sols,snpNameList);
    outfile(genos_file,snpNameList, elite_sols);

    return 0;
}

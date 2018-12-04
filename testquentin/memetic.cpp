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

#include "CSVParser.hpp"
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

struct parents_pairs {
  patternscore parent1;
  patternscore parent2;
};

vector<string> get_snp_list(string genos_file){
  vector<string> tokens;
  ifstream file;
  file.open(genos_file);
  string SNPheader;
  getline(file,SNPheader);
  file.close();
  boost::algorithm::split(tokens, SNPheader, boost::is_any_of(","));
  return tokens;
}

float g_test_2SNP(contingence2SNP contingence2){
  int const nbrcolonnes(9);
  int const nbrligne(2);
  //cout<<nbrcolonnes<<endl;

  for(int j=0; j<(nbrcolonnes); ++j)
  {    for(int i=0; i<(nbrligne); ++i)
      {

         //cout<<contingence2[nbrligne][j]<<"+"<<contingence2[i][j]<<"=";
         contingence2[nbrligne][j] += contingence2[i][j];
         //cout<<contingence2[nbrligne][j]<<" ";

      }
      //cout<<";"<<endl;
  }

  for(int j=0; j<(nbrcolonnes); ++j)
  {    for(int i=0; i<(nbrligne); ++i)
      {
         contingence2[i][nbrcolonnes] += contingence2[i][j];
      }
  }
  for(int i(0); i<(nbrligne); ++i)
      {
         contingence2[nbrligne][nbrcolonnes] += contingence2[i][nbrcolonnes];
      }

  for (int i=0;i<3;i++){
    for (int j=0;j<10;j++){
      cout<<contingence2[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;

  float contingencetheo[nbrligne][nbrcolonnes];
  for (int i=0;i<(nbrligne);i++){ //Set contingency matrix to 0
    for (int j=0;j<(nbrcolonnes);j++){
      contingencetheo[i][j]=0;
    }
  }
  for(int i=0; i<(nbrligne); ++i){
    for(int j=0; j<(nbrcolonnes); ++j){
       //cout<<contingence2[i][nbrcolonnes]<<"*"<<contingence2[nbrligne][j]<<"/"<<contingence2[nbrligne][nbrcolonnes];
      contingencetheo[i][j] = (float) contingence2[i][nbrcolonnes]* (float)contingence2[nbrligne][j]/(float)contingence2[nbrligne][nbrcolonnes];
      // cout<<"="<<contingencetheo[i][j]<<endl;
    }
  }
  for (int i=0;i<2;i++){
    for (int j=0;j<9;j++){
      cout<<contingencetheo[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;

float scorekhi2=0;
/*
  for(int i(0); i<(nbrligne); ++i)
  {    for(int j(0); j<(nbrcolonnes); ++j)
      {
         scorekhi2 += (pow(((contingence2[i][j]-(contingencetheo[i][j]))),2)/contingencetheo[i][j]);
         //cout<<scorekhi2<<"="<<contingence2[i][j]<<"-"<<contingencetheo[i][j]<<"^2"<<"/"<<contingencetheo[i][j]<<endl;
      }
  }
  */

int test;
test=0;
/*
unsigned ncells = nbrligne*nbrcolonnes;
int count_inf_5 = 0;
for(unsigned i=0; i<nbrligne; ++i){
  for(unsigned j=0; j<nbrcolonnes; ++j){
    if(contingencetheo[i][j] < 5){
      count_inf_5 ++;
      if((double)count_inf_5 / ncells > 0.2){
            test=1;
      }
    }
  }
}
for(unsigned i=0; i<nbrligne; ++i){
  for(unsigned j=0; j<nbrcolonnes; ++j){
    if(contingence2[i][j] < 5){
      test=1;
    }
  }
}
if(test==1){
  cout<<"The test isn't reliable"<<endl;
}
*/
  int df;
  df=(nbrligne-1)*(nbrcolonnes-1);
  float pval;
  if(test==0){
    for(int i(0); i<(nbrligne); ++i){
      for(int j(0); j<(nbrcolonnes); ++j){
        if (contingence2[i][j] != 0 ){
          double div=(double)contingence2[i][j]/contingencetheo[i][j];
          scorekhi2 += contingence2[i][j] * log(div);
        }
      }
    }
    scorekhi2  *= 2;
    boost::math::chi_squared_distribution<double> chi2_dist(df);
    pval = 1 - boost::math::cdf(chi2_dist, scorekhi2);
    if(pval == 0){
        pval = 2.0e-16;}
    cout<<"score: "<<scorekhi2<<endl;
    cout<<"p: "<<pval<<endl;
  }
  else{
    scorekhi2=0;
    pval=1;
    cout<<"score: "<<scorekhi2<<endl;
    cout<<"p: "<<pval<<endl;
  }
  return(scorekhi2);
}

float g_test_3SNP(contingence3SNP contingence2){
  int const nbrcolonnes(27);
  int const nbrligne(2);
  //cout<<nbrcolonnes<<endl;

  for(int j=0; j<(nbrcolonnes); ++j)
  {    for(int i=0; i<(nbrligne); ++i)
      {

         //cout<<contingence2[nbrligne][j]<<"+"<<contingence2[i][j]<<"=";
         contingence2[nbrligne][j] += contingence2[i][j];
         //cout<<contingence2[nbrligne][j]<<" ";

      }
      //cout<<";"<<endl;
  }

  for(int j=0; j<(nbrcolonnes); ++j)
  {    for(int i=0; i<(nbrligne); ++i)
      {
         contingence2[i][nbrcolonnes] += contingence2[i][j];
      }
  }
  for(int i(0); i<(nbrligne); ++i)
      {
         contingence2[nbrligne][nbrcolonnes] += contingence2[i][nbrcolonnes];
      }

  for (int i=0;i<3;i++){
    for (int j=0;j<28;j++){
      cout<<contingence2[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;

  float contingencetheo[nbrligne][nbrcolonnes];
  for (int i=0;i<(nbrligne);i++){ //Set contingency matrix to 0
    for (int j=0;j<(nbrcolonnes);j++){
      contingencetheo[i][j]=0;
    }
  }
  for(int i=0; i<(nbrligne); ++i){
    for(int j=0; j<(nbrcolonnes); ++j){
       //cout<<contingence2[i][nbrcolonnes]<<"*"<<contingence2[nbrligne][j]<<"/"<<contingence2[nbrligne][nbrcolonnes];
      contingencetheo[i][j] = (float) contingence2[i][nbrcolonnes]* (float)contingence2[nbrligne][j]/(float)contingence2[nbrligne][nbrcolonnes];
      // cout<<"="<<contingencetheo[i][j]<<endl;
    }
  }
  for (int i=0;i<2;i++){
    for (int j=0;j<27;j++){
      cout<<contingencetheo[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;

float scorekhi2=0;
/*
  for(int i(0); i<(nbrligne); ++i)
  {    for(int j(0); j<(nbrcolonnes); ++j)
      {
         scorekhi2 += (pow(((contingence2[i][j]-(contingencetheo[i][j]))),2)/contingencetheo[i][j]);
         //cout<<scorekhi2<<"="<<contingence2[i][j]<<"-"<<contingencetheo[i][j]<<"^2"<<"/"<<contingencetheo[i][j]<<endl;
      }
  }
  */

int test;
test=0;
/*
unsigned ncells = nbrligne*nbrcolonnes;
int count_inf_5 = 0;
for(unsigned i=0; i<nbrligne; ++i){
  for(unsigned j=0; j<nbrcolonnes; ++j){
    if(contingencetheo[i][j] < 5){
      count_inf_5 ++;
      if((double)count_inf_5 / ncells > 0.2){
            test=1;
      }
    }
  }
}
for(unsigned i=0; i<nbrligne; ++i){
  for(unsigned j=0; j<nbrcolonnes; ++j){
    if(contingence2[i][j] < 5){
      test=1;
    }
  }
}
if(test==1){
  cout<<"The test isn't reliable"<<endl;
}
*/
  int df;
  df=(nbrligne-1)*(nbrcolonnes-1);
  float pval;
  if(test==0){
    for(int i(0); i<(nbrligne); ++i){
      for(int j(0); j<(nbrcolonnes); ++j){
        if (contingence2[i][j] != 0 ){
          double div=(double)contingence2[i][j]/contingencetheo[i][j];
          scorekhi2 += contingence2[i][j] * log(div);
        }
      }
    }
    scorekhi2  *= 2;
    boost::math::chi_squared_distribution<double> chi2_dist(df);
    pval = 1 - boost::math::cdf(chi2_dist, scorekhi2);
    if(pval == 0){
        pval = 2.0e-16;}
    cout<<"score: "<<scorekhi2<<endl;
    cout<<"p: "<<pval<<endl;
  }
  else{
    scorekhi2=0;
    pval=1;
    cout<<"score: "<<scorekhi2<<endl;
    cout<<"p: "<<pval<<endl;
  }
  return(scorekhi2);
}






void create_contingency_table_pattern2(int l1,int l2,contingence2SNP* adr_contingence,blas_matrix genos, blas_matrix phenos_m){
  int k;
  for (int i=0;i<3;i++){ //Set contingency matrix to 0
    for (int j=0;j<10;j++){
      (*adr_contingence)[i][j]=0;
    }
  }
  for (k=0;k<int(genos.size1());k++){ //For each individuals in the data imported
    //Filling the contingency matrix
    if (phenos_m(k,0)==1){
      if (genos(k,l1)==0){
        if (genos(k,l2)==0){
          (*adr_contingence)[0][0]+=1;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[0][1]+=1;
        }
        else{
          (*adr_contingence)[0][2]+=1;
        }
      }
      else if (genos(k,l1)==1){
        if (genos(k,l2)==0){
          (*adr_contingence)[0][3]+=1;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[0][4]+=1;
        }
        else{
          (*adr_contingence)[0][5]+=1;
        }
      }
      else{
        if (genos(k,l2)==0){
          (*adr_contingence)[0][6]+=1;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[0][7]+=1;
        }
        else{
          (*adr_contingence)[0][8]+=1;
        }
      }
    }
    else{
      if (genos(k,l1)==0){
        if (genos(k,l2)==0){
          (*adr_contingence)[1][0]+=1;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[1][1]+=1;
        }
        else{
          (*adr_contingence)[1][2]+=1;
        }
      }
      else if (genos(k,l1)==1){
        if (genos(k,l2)==0){
          (*adr_contingence)[1][3]+=1;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[1][4]+=1;
        }
        else{
          (*adr_contingence)[1][5]+=1;
        }
      }
      else{
        if (genos(k,l2)==0){
          (*adr_contingence)[1][6]+=1;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[1][7]+=1;
        }
        else{
          (*adr_contingence)[1][8]+=1;
        }
      }
    }
  }
}

void create_contingency_table_pattern3(int l1,int l2,int l3,contingence3SNP* adr_contingence,blas_matrix genos, blas_matrix phenos_m){
  int k;
  for (int i=0;i<3;i++){ //Set contingency matrix to 0
    for (int j=0;j<28;j++){
      (*adr_contingence)[i][j]=0;
    }
  }
  for (k=0;k<int(genos.size1());k++){ //For each individuals in the data imported
    //Filling the contingency matrix
    if (phenos_m(k,0)==1){
      if (genos(k,l1)==0){
        if (genos(k,l2)==0){
          if (genos(k,l3)==0){
            (*adr_contingence)[0][0]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][1]+=1;
          }
          else{
            (*adr_contingence)[0][2]+=1;
          }
        }
        else if (genos(k,l2)==1){
          if (genos(k,l3)==0){
            (*adr_contingence)[0][3]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][4]+=1;
          }
          else{
            (*adr_contingence)[0][5]+=1;
          }
        }
        else{
          if (genos(k,l3)==0){
            (*adr_contingence)[0][6]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][7]+=1;
          }
          else{
            (*adr_contingence)[0][8]+=1;
          }
        }
      }
      else if (genos(k,l1)==1){
        if (genos(k,l2)==0){
          if (genos(k,l3)==0){
            (*adr_contingence)[0][9]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][10]+=1;
          }
          else{
            (*adr_contingence)[0][11]+=1;
          }
        }
        else if (genos(k,l2)==1){
          if (genos(k,l3)==0){
            (*adr_contingence)[0][12]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][13]+=1;
          }
          else{
            (*adr_contingence)[0][14]+=1;
          }
        }
        else{
          if (genos(k,l3)==0){
            (*adr_contingence)[0][15]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][16]+=1;
          }
          else{
            (*adr_contingence)[0][17]+=1;
          }
        }
      }
      else {
        if (genos(k,l2)==0){
          if (genos(k,l3)==0){
            (*adr_contingence)[0][18]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][19]+=1;
          }
          else{
            (*adr_contingence)[0][20]+=1;
          }
        }
        else if (genos(k,l2)==1){
          if (genos(k,l3)==0){
            (*adr_contingence)[0][21]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][22]+=1;
          }
          else{
            (*adr_contingence)[0][23]+=1;
          }
        }
        else{
          if (genos(k,l3)==0){
            (*adr_contingence)[0][24]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[0][25]+=1;
          }
          else{
            (*adr_contingence)[0][26]+=1;
          }
        }
      }
    }
    else{
      if (genos(k,l1)==0){
        if (genos(k,l2)==0){
          if (genos(k,l3)==0){
            (*adr_contingence)[1][0]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][1]+=1;
          }
          else{
            (*adr_contingence)[1][2]+=1;
          }
        }
        else if (genos(k,l2)==1){
          if (genos(k,l3)==0){
            (*adr_contingence)[1][3]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][4]+=1;
          }
          else{
            (*adr_contingence)[1][5]+=1;
          }
        }
        else{
          if (genos(k,l3)==0){
            (*adr_contingence)[1][6]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][7]+=1;
          }
          else{
            (*adr_contingence)[1][8]+=1;
          }
        }
      }
      else if (genos(k,l1)==1){
        if (genos(k,l2)==0){
          if (genos(k,l3)==0){
            (*adr_contingence)[1][9]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][10]+=1;
          }
          else{
            (*adr_contingence)[1][11]+=1;
          }
        }
        else if (genos(k,l2)==1){
          if (genos(k,l3)==0){
            (*adr_contingence)[1][12]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][13]+=1;
          }
          else{
            (*adr_contingence)[1][14]+=1;
          }
        }
        else{
          if (genos(k,l3)==0){
            (*adr_contingence)[1][15]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][16]+=1;
          }
          else{
            (*adr_contingence)[1][17]+=1;
          }
        }
      }
      else {
        if (genos(k,l2)==0){
          if (genos(k,l3)==0){
            (*adr_contingence)[1][18]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][19]+=1;
          }
          else{
            (*adr_contingence)[1][20]+=1;
          }
        }
        else if (genos(k,l2)==1){
          if (genos(k,l3)==0){
            (*adr_contingence)[1][21]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][22]+=1;
          }
          else{
            (*adr_contingence)[1][23]+=1;
          }
        }
        else{
          if (genos(k,l3)==0){
            (*adr_contingence)[1][24]+=1;
          }
          else if (genos(k,l3)==1){
            (*adr_contingence)[1][25]+=1;
          }
          else{
            (*adr_contingence)[1][26]+=1;
          }
        }
      }
    }
  }
}

vector<parents_pairs> select_pairs_of_individuals_to_be_crossed(vector<patternscore> n_pairs_selected_parents){
  parents_pairs couple;
  vector<parents_pairs>pairs_of_parents;
  vector<patternscore>nodes(n_pairs_selected_parents.begin(),n_pairs_selected_parents.end());
  srand(time(0));
  random_shuffle(nodes.begin(), nodes.end());
  for (int i=0;i<n_pairs_selected_parents.size();i+=2){
    if(i+2>n_pairs_selected_parents.size() || i+1>n_pairs_selected_parents.size()){
      break;
    }
    couple.parent1=nodes[i];
    couple.parent2=nodes[i+1];
    pairs_of_parents.push_back(couple);
  }return(pairs_of_parents);
}

vector<patternscore> create_two_children_for_each_selected_pair_of_parents(vector<parents_pairs> pairs_of_parents){
   patternscore son;
   patternscore daugther;
   vector<patternscore> children_parents;
   for (int i=0;i<pairs_of_parents.size();++i){
     if(pairs_of_parents[i].parent1.pattern1==pairs_of_parents[i].parent2.pattern1){continue;}
     son.pattern1=pairs_of_parents[i].parent1.pattern1;
     son.pattern2=pairs_of_parents[i].parent2.pattern1;
     children_parents.push_back(son);
     if(pairs_of_parents[i].parent1.pattern2==pairs_of_parents[i].parent2.pattern2){continue;}
     daugther.pattern1=pairs_of_parents[i].parent1.pattern2;
     daugther.pattern2=pairs_of_parents[i].parent2.pattern2;
     children_parents.push_back(daugther);
   } return(children_parents);

 }

void perform_one_mutation_per_child(vector<patternscore>children_parents,int prob_mutation){
  for (int i=0;i<children_parents.size();++i){
    int mutation = rand() % 100;
    if (mutation < prob_mutation){
      break;
    }else{
      int type = rand() % 2 + 0;
      int parentpattern = rand() % 1 + 0;
      int mutpattern = rand() % 1 + 0;
      switch (type) {
        case 0:
          break;
        case 1:
          {int snp = rand() % children_parents.size();
          if(parentpattern ==0 ){if(mutpattern ==0){
            children_parents[i].pattern1=children_parents[snp].pattern1;}
            else{children_parents[i].pattern1=children_parents[snp].pattern2;}}
          else{if(mutpattern ==0){
            children_parents[i].pattern2=children_parents[snp].pattern1;}
            else{children_parents[i].pattern2=children_parents[snp].pattern2;}}}break;
        case 2:
          break;
        }
      }
    }
  }

  void cout_list_indiv(vector<parents_pairs> list_to_cout_indiv){
    for (vector<parents_pairs>::iterator it=list_to_cout_indiv.begin();it!=list_to_cout_indiv.end();it++){
        cout<<"Parent"<<endl;
        cout<<"parent1:"<<(*it).parent1.pattern1<<","<<(*it).parent1.pattern2<<"   "<<(*it).parent1.score<<endl;
        cout<<"parent2:"<<(*it).parent2.pattern1<<","<<(*it).parent2.pattern2<<"   "<<(*it).parent2.score<<endl;
      }
    }

    void cout_list(vector<patternscore> list_to_cout){
      for (vector<patternscore>::iterator it=list_to_cout.begin();it!=list_to_cout.end();it++){
        if ((*it).pattern3==""){
          cout<<(*it).pattern1<<","<<(*it).pattern2<<endl;
          cout<<(*it).score<<endl;
        }
        else{
          cout<<(*it).pattern1<<","<<(*it).pattern2<<(*it).pattern3<<endl;
          cout<<(*it).score<<endl;
        }
      }
    }

    int calculate_delta(patternscore s, patternscore sB){
      int diff=0;
      if (s.pattern1!=sB.pattern1 && s.pattern1!=sB.pattern2 && s.pattern1!=sB.pattern3){
        diff=diff+1;
      }
      if (s.pattern2!=sB.pattern1 && s.pattern2!=sB.pattern2 && s.pattern2!=sB.pattern3){
        diff=diff+1;
      }
      if (s.pattern3!=sB.pattern1 && s.pattern3!=sB.pattern2 && s.pattern3!=sB.pattern3){
        diff=diff+1;
      }
      return(diff);
    }



    vector<patternscore> neighbours(patternscore s,vector<patternscore> patternscoreList){
      vector<patternscore> s_neighbours;
      for (int i=0;i<patternscoreList.size();i++){
        int delta=calculate_delta(s,patternscoreList[i]);
        if (delta==1){
          s_neighbours.push_back(patternscoreList[i]);
        }
      }
      return(s_neighbours);
    }

    patternscore select_closest_neighbor_to_guiding_solution(patternscore s,patternscore sB, vector<patternscore> patternscoreList){
      patternscore s_closest_neighbour;
      vector<patternscore> s_neighbours=neighbours(s,patternscoreList);
      int min_delta=999;
      for (int i=0;i<s_neighbours.size();i++){
        int delta=calculate_delta(s_neighbours[i],sB);
        if (delta<min_delta){
          min_delta=delta;
          s_closest_neighbour=s_neighbours[i];
        }
      }
      return(s_closest_neighbour);
    }


    patternscore hill_climbing_lc(patternscore s_closest_neighbour, vector<patternscore> patternscoreList){
      vector<patternscore> s_neighbours = neighbours(s_closest_neighbour,patternscoreList);
      patternscore actual_s = s_closest_neighbour;
      for (int i=0;i<patternscoreList.size();i++){
        if (s_neighbours[i].score>actual_s.score){
          actual_s=s_neighbours[i];
        }
      }
      return(actual_s);
    }

    vector<patternscore> initialize_population(int n,vector<patternscore> patternscoreList){
      vector<patternscore> pop;
      vector<patternscore> nodes(patternscoreList.begin(),patternscoreList.end());
      srand(time(0));
      random_shuffle(nodes.begin(), nodes.end());
      for (vector<patternscore>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
        if(pop.size()<n){
          pop.push_back(*it);
        }
        else{
          break;
        }
      }
      return(pop);
    }


    float add_gtest_score (patternscore pattern,blas_matrix genos,blas_matrix phenos_m){
      float score;
      if(pattern.snp3==NULL){
        contingence2SNP contingence2;
        contingence2SNP* adr_contingence2 = &contingence2;
        create_contingency_table_pattern2(pattern.snp1,pattern.snp2,adr_contingence2,genos,phenos_m);
        score=g_test_2SNP(contingence2);
      }
      else{
        contingence3SNP contingence3;
        contingence3SNP* adr_contingence3 = &contingence3;
        create_contingency_table_pattern3(pattern.snp1,pattern.snp2,pattern.snp3,adr_contingence3,genos,phenos_m);
        score=g_test_3SNP(contingence3);
      }
      return(score);
    }

    int main()
    {
        // Arguments
        string genos_file;
        genos_file="simu2_Genotype_1.csv";
        string phenos_file;
        phenos_file="simu2_Phenotype_1.csv";
        int header = 1;
        char separator = ',';

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
        /*
                create_contingency_table_pattern2(l1,l2,adr_contingence2,genos,phenos_m);
                cout<<snpNameList[l1]<<","<<snpNameList[l2]<<endl;
                int countNonStat=0;
                for (int i=0;i<2;i++){
                  for (int j=0;j<9;j++){
                    cout<<contingence2[i][j]<<" ";
                    if (contingence2[i][j]<5){
                      countNonStat++;
                    }
                  }
                  cout<<endl;
                }
                cout<<countNonStat<<" valeurs dans le tableau inférieures à 5."<<endl;
                cout<<endl;
                float scorekhi2=g_test_2SNP(contingence2);
        */

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
        /*
                  create_contingency_table_pattern3(l1,l2,l3,adr_contingence3,genos,phenos_m);
                  cout<<snpNameList[l1]<<","<<snpNameList[l2]<<","<<snpNameList[l3]<<endl;
                  int countNonStat=0;
                  for (int i=0;i<2;i++){
                    for (int j=0;j<27;j++){
                      cout<<contingence3[i][j]<<" ";
                      if (contingence3[i][j]<5){
                        countNonStat++;
                      }
                    }
                    cout<<endl;
                  }
                  cout<<countNonStat<<" valeurs dans le tableau inférieures à 5."<<endl;
                  cout<<endl;

                  float scorekhi2=g_test_3SNP(contingence3);
        */
                }
              }
            }





        int prob_mutation = 1;
        int n_it = 3;
        int n = 10;
        int j = 1;
        vector<patternscore>n_pairs_selected_parents=initialize_population(n, patternscoreList);


        for (int i=0;i<n_pairs_selected_parents.size();i++){
            n_pairs_selected_parents[i].score=add_gtest_score(n_pairs_selected_parents[i],genos,phenos_m);
        }

        cout<<"lol"<<endl;
        while (j < n_it){
          cout<<"darons:"<<endl;
          vector<parents_pairs> pairs_of_parents=select_pairs_of_individuals_to_be_crossed(n_pairs_selected_parents);
          cout_list_indiv(pairs_of_parents);
          cout<<"enfant:"<<endl;
          vector<patternscore> children_parents;
          children_parents=create_two_children_for_each_selected_pair_of_parents(pairs_of_parents);
          cout_list(children_parents);


          for (int i=0;i<children_parents.size();i++){
              cout<<children_parents.size()<<"i"<<endl;
              float score=add_gtest_score(children_parents[1],genos,phenos_m);
          }
          cout_list(children_parents);
          perform_one_mutation_per_child(children_parents,prob_mutation);

          for (int i=0;i<children_parents.size();i++){
              children_parents[i].score=add_gtest_score(children_parents[i],genos,phenos_m);
          }
          j=j+1;
        }
        return 0;
    }

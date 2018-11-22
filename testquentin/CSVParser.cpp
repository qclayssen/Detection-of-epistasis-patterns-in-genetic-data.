//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
//Parsing of CSV Files for Epistasis detection
//Created :
//Modified :

#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <boost/math/distributions/chi_squared.hpp>

#include "CSVParser.hpp"
#include "common.h"

using namespace std;
using namespace std::chrono;

typedef int contingence2SNP[3][10];

void create_contingency_table(int nb_sample,contingence2SNP* adr_contingence,blas_matrix genos, blas_matrix phenos_m){
  int k=nb_sample;
  int l1,l2;
  for (int i=0;i<3;i++){ //Set contingency matrix to 0
    for (int j=0;j<10;j++){
      (*adr_contingence)[i][j]=0;
    }
  }
  for (l1=0;l1<int(genos.size1())-2;l1++){ //First SNP of the pattern
    for (l2=1;l2<int(genos.size1())-1;l2++){ //Second SNP of the pattern
      //Filling the contingency matrix
      if (phenos_m(l1,0)==1){
        if (genos(l1,k)==0){
          if (genos(l2,k)==0){
            (*adr_contingence)[0][0]+=1;
          }
          else if (genos(l2,k)==1){
            (*adr_contingence)[0][1]+=1;
          }
          else{
            (*adr_contingence)[0][2]+=1;
          }
        }
        else if (genos(l1,k)==1){
          if (genos(l2,k)==0){
            (*adr_contingence)[0][3]+=1;
          }
          else if (genos(l2,k)==1){
            (*adr_contingence)[0][4]+=1;
          }
          else{
            (*adr_contingence)[0][5]+=1;
          }
        }
        else{
          if (genos(l2,k)==0){
            (*adr_contingence)[0][6]+=1;
          }
          else if (genos(l2,k)==1){
            (*adr_contingence)[0][7]+=1;
          }
          else{
            (*adr_contingence)[0][8]+=1;
          }
        }
      }
      else{
        if (genos(l1,k)==0){
          if (genos(l2,k)==0){
            (*adr_contingence)[1][0]++;
          }
          else if (genos(l2,k)==1){
            (*adr_contingence)[1][1]++;
          }
          else{
            (*adr_contingence)[1][2]++;
          }
        }
        else if (genos(l1,k)==1){
          if (genos(l2,k)==0){
            (*adr_contingence)[1][3]++;
          }
          else if (genos(l2,k)==1){
            (*adr_contingence)[1][4]++;
          }
          else{
            (*adr_contingence)[1][5]++;
          }
        }
        else{
          if (genos(l2,k)==0){
            (*adr_contingence)[1][6]++;
          }
          else if (genos(l2,k)==1){
            (*adr_contingence)[1][7]++;
          }
          else{
            (*adr_contingence)[1][8]++;
          }
        }
      }
    }
  }
}

int main()
{
    // Arguments
    string genos_file;
    genos_file="simu1_Genotype_1.csv";
    cout<<"Genotype file path : "<<endl;
    string phenos_file;
    phenos_file="simu1_Phenotype_1.csv";
    cout<<"Phenotype file path : "<<endl;
    int header = 1;
    char separator = ',';

//  DATA IMPORTATION
    CSVParser<int> genos_csv(genos_file, separator, header);
    CSVParser<int> phenos_csv(phenos_file, separator, header);
    blas_matrix genos = genos_csv.data();
    blas_matrix phenos_m = phenos_csv.data();
    blas_column phenos(phenos_m, 0);
    cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;
    int k;

    contingence2SNP contingence;
    contingence2SNP* adr_contingence = &contingence;
    for (k=0;k<int(genos.size2())-1;k++){ //For each individuals in the data imported
      create_contingency_table(k,adr_contingence,genos,phenos_m);
      int countNonStat=0;
      for (int i=0;i<2;i++){
        for (int j=0;j<9;j++){
          cout<<contingence[i][j]<<" ";
          if (contingence[i][j]<5){
            countNonStat++;
          }
        }
        cout<<endl;
      }
      cout<<endl;
      cout<<countNonStat<<" valeurs dans le tableau inférieures à 5."<<endl;

      int const nbrcolonnes(9);
      int const nbrligne(2);
      //cout<<nbrcolonnes<<endl;

      for(int j(0); j<(nbrcolonnes); ++j)
      {    for(int i(0); i<(nbrligne); ++i)
          {

             //cout<<contingence[nbrligne][j]<<"+"<<contingence[i][j]<<"=";
             contingence[nbrligne][j] += contingence[i][j];
             //cout<<contingence[nbrligne][j]<<" ";

          }
          //cout<<";"<<endl;
      }
      for (int i=0;i<3;i++){
        for (int j=0;j<10;j++){
          cout<<contingence[i][j]<<" ";
        }
        cout<<endl;
      }
      cout<<endl;



      for(int j(0); j<(nbrcolonnes); ++j)
      {    for(int i(0); i<(nbrligne); ++i)
          {
             contingence[i][nbrcolonnes] += contingence[i][j];
          }
      }
      for(int i(0); i<(nbrligne); ++i)
          {
             contingence[nbrligne][nbrcolonnes] += contingence[i][nbrcolonnes];
          }
      for (int i=0;i<3;i++){
        for (int j=0;j<10;j++){
          cout<<contingence[i][j]<<" ";
        }
        cout<<endl;
      }
      cout<<endl;



      int contingencetheo[nbrcolonnes][nbrligne];
      for(int i(0); i<(nbrligne); ++i)
      {    for(int j(0); j<(nbrcolonnes); ++j)
          {
             //cout<<contingence[i][nbrcolonnes]<<"*"<<contingence[nbrligne][j]<<"/"<<contingence[nbrligne][nbrcolonnes];
             contingencetheo[i][j] = contingence[i][nbrcolonnes]*contingence[nbrligne][j]/contingence[nbrligne][nbrcolonnes];
            // cout<<"="<<contingencetheo[i][j]<<endl;
          }
      }    for (int i=0;i<2;i++){
            for (int j=0;j<9;j++){
              cout<<contingencetheo[i][j]<<" ";
            }
            cout<<endl;
          }
          cout<<endl;


          float scorekhi2=0;
          for(int i(0); i<(nbrligne); ++i)
          {    for(int j(0); j<(nbrcolonnes); ++j)
              {
                      double div = (double) contingencetheo[i][j] / contingence[i][j];
                      scorekhi2 += contingencetheo[i][j] * log(div);

              }

          }
          scorekhi2  *= 2;


          float pval;
          cout<<"score: "<<scorekhi2<<endl;
          boost::math::chi_squared mydist(8);
          pval = 1 - boost::math::cdf(mydist, scorekhi2);
          cout<<"p: "<<pval<<endl;
          return(scorekhi2);



        }
    return 0;
  }

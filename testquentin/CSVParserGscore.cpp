//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
//Parsing of CSV Files for Epistasis detection
//Created :
//Modified :

#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/algorithm/string.hpp>
#include <list>
#include "CSVParser.hpp"
#include "common.h"

using namespace std;
using namespace std::chrono;

typedef int contingence2SNP[3][10];
struct patternscore {
  string pattern1;
  string pattern2;
  float score;
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

void create_contingency_table(int l1,int l2,contingence2SNP* adr_contingence,blas_matrix genos, blas_matrix phenos_m){
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
          (*adr_contingence)[1][0]++;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[1][1]++;
        }
        else{
          (*adr_contingence)[1][2]++;
        }
      }
      else if (genos(k,l1)==1){
        if (genos(k,l2)==0){
          (*adr_contingence)[1][3]++;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[1][4]++;
        }
        else{
          (*adr_contingence)[1][5]++;
        }
      }
      else{
        if (genos(k,l2)==0){
          (*adr_contingence)[1][6]++;
        }
        else if (genos(k,l2)==1){
          (*adr_contingence)[1][7]++;
        }
        else{
          (*adr_contingence)[1][8]++;
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
    string phenos_file;
    phenos_file="simu1_Phenotype_1.csv";
    int header = 1;
    char separator = ',';

//  DATA IMPORTATION
    CSVParser<int> genos_csv(genos_file, separator, header);
    CSVParser<int> phenos_csv(phenos_file, separator, header);
    blas_matrix genos = genos_csv.data();
    blas_matrix phenos_m = phenos_csv.data();
    blas_column phenos(phenos_m, 0);
    cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;
    int l1,l2;

    vector<string> snpNameList;
    snpNameList = get_snp_list(genos_file);

    list<patternscore> patternscoreList;
    list<patternscore>::iterator iterpatternscoreList;

    contingence2SNP contingence;
    contingence2SNP* adr_contingence = &contingence;
    for (l1=0;l1<int(genos.size2())-2;l1++){ //First SNP of the pattern
      for (l2=l1+1;l2<int(genos.size2())-1;l2++){ //Second SNP of the pattern
        create_contingency_table(l1,l2,adr_contingence,genos,phenos_m);
        cout<<snpNameList[l1]<<","<<snpNameList[l2]<<endl;
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
        cout<<countNonStat<<" valeurs dans le tableau inférieures à 5."<<endl;
        cout<<endl;

        int const nbrcolonnes(9);
        int const nbrligne(2);
        //cout<<nbrcolonnes<<endl;

        for(int j=0; j<(nbrcolonnes); ++j)
        {    for(int i=0; i<(nbrligne); ++i)
            {

               //cout<<contingence[nbrligne][j]<<"+"<<contingence[i][j]<<"=";
               contingence[nbrligne][j] += contingence[i][j];
               //cout<<contingence[nbrligne][j]<<" ";

            }
            //cout<<";"<<endl;
        }

        for(int j=0; j<(nbrcolonnes); ++j)
        {    for(int i=0; i<(nbrligne); ++i)
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
        for(int i=0; i<(nbrligne); ++i)
        {    for(int j=0; j<(nbrcolonnes); ++j)
            {
               //cout<<contingence[i][nbrcolonnes]<<"*"<<contingence[nbrligne][j]<<"/"<<contingence[nbrligne][nbrcolonnes];
               contingencetheo[i][j] = contingence[i][nbrcolonnes]*contingence[nbrligne][j]/contingence[nbrligne][nbrcolonnes];
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



              float scorekhi2=0;/*
                for(int i(0); i<(nbrligne); ++i)
                {    for(int j(0); j<(nbrcolonnes); ++j)
                    {
                       scorekhi2 += (pow(((contingence[i][j]-(contingencetheo[i][j]))),2)/contingencetheo[i][j]);
                       //cout<<scorekhi2<<"="<<contingence[i][j]<<"-"<<contingencetheo[i][j]<<"^2"<<"/"<<contingencetheo[i][j]<<endl;
                    }
                }*/
              int test;
              test=0;
              unsigned ncells = nbrligne*nbrcolonnes;
              int count_inf_5 = 0;
              for(unsigned i=0; i<nbrligne; ++i)
              {
                  for(unsigned j=0; j<nbrcolonnes; ++j)
                  {
                      if(contingencetheo[i][j] < 0 || contingencetheo[i][j]!=contingencetheo[i][j]) // test for nan
                      {
                          test=1;
                      }
                      if(contingencetheo[i][j] < 5)
                      {
                          count_inf_5 ++;
                          if((double)count_inf_5 / ncells > 0.2)
                          {
              //                    cout << "Not reliable test, ";
              //                    cout << "expected: " << e << endl;
                              test=1;
                          }
                      }
                  }
                }
                for(unsigned i=0; i<nbrligne; ++i)
                {
                   for(unsigned j=0; j<nbrcolonnes; ++j)
                   {
                         if(contingence[i][j] < 5)
                         test=1;
                   }
                }
                cout<<"test:"<<test<<endl;



                float pval;
                if(test==0)
                {
                for(int i(0); i<(nbrligne); ++i)
                {    for(int j(0); j<(nbrcolonnes); ++j)
                    {
                      if (contingence[i][j] != 0 ){
                            double div = (double) contingence[i][j] / contingencetheo[i][j];
                            cout<<"contingence[i][j]:"<<contingence[i][j]<<endl;
                            cout<<"contingencetheo[i][j]:"<<contingencetheo[i][j]<<endl;
                            cout<<"div:"<<div<<endl;
                            cout<<"log div:"<<log(div)<<endl;
                            cout<<"scorekhi2:"<<scorekhi2<<endl;
                            scorekhi2 += contingence[i][j] * log(div);
                          }

                    }
                  }            scorekhi2  *= 2;
                              boost::math::chi_squared mydist(8);
                              pval = 1 - boost::math::cdf(mydist, scorekhi2);
                              if(pval == 0){
                                  pval = 2.0e-16;}
                              cout<<"score: "<<scorekhi2<<endl;
                              cout<<"p: "<<pval<<endl;



                }
                else {
                  scorekhi2=0;
                  pval=1;
                  cout<<"score: "<<endl;
                  cout<<"p: "<<endl;
                  continue;         }
                  patternscore p1;
                  p1.pattern1=snpNameList[l1];
                  p1.pattern2=snpNameList[l2];
                  p1.score=scorekhi2;
                  patternscoreList.push_back(p1);

                }
                  }for (iterpatternscoreList=patternscoreList.begin();iterpatternscoreList!=patternscoreList.end();iterpatternscoreList++){cout<<(*iterpatternscoreList).pattern1<<","<<(*iterpatternscoreList).pattern2<<"="<<(*iterpatternscoreList).score<<endl;}

                          return 0;
        }

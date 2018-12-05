#include <iostream>
#include <string>
#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;
namespace blas = boost::numeric::ublas;

struct patternscore {
  string pattern1;
  string pattern2;
  float score;
};

patternscore gscore(){
  for (l1=0;l1<int(genos.size2())-1;l1++){ //First SNP of the pattern
    for (l2=l1+1;l2<int(genos.size2());l2++){ //Second SNP of the pattern
      create_contingency_table_pattern2(l1,l2,adr_contingence2,genos,phenos_m);
      //cout<<snpNameList[l1]<<","<<snpNameList[l2]<<endl;
      int countNonStat=0;
      for (int i=0;i<2;i++){
        for (int j=0;j<9;j++){
          cout<<contingence2[i][j]<<";";
          if (contingence2[i][j]<5){
            countNonStat++;
          }
        }
        cout<<endl;
      }
      //cout<<countNonStat<<" valeurs dans le tableau inférieures à 5."<<endl;
      cout<<endl;

      int const nbrcolonnes(9);
      int const nbrligne(2);
      ////cout<<nbrcolonnes<<endl;
      for(int j=0; j<(nbrcolonnes); ++j)
      {    for(int i=0; i<(nbrligne); ++i)
          {
             ////cout<<contingence2[nbrligne][j]<<"+"<<contingence2[i][j]<<"=";
             contingence2[nbrligne][j] += contingence2[i][j];
             ////cout<<contingence2[nbrligne][j]<<" ";
          }
          ////cout<<";"<<endl;
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
          //cout<<contingence2[i][j]<<" ";
        }
        //cout<<endl;
      }
      //cout<<endl;
      float contingencetheo[nbrligne][nbrcolonnes];
      for(int i=0; i<(nbrligne); ++i)
      {    for(int j=0; j<(nbrcolonnes); ++j)
          {
             ////cout<<contingence2[i][nbrcolonnes]<<"*"<<contingence2[nbrligne][j]<<"/"<<contingence2[nbrligne][nbrcolonnes];
             contingencetheo[i][j] = contingence2[i][nbrcolonnes]*contingence2[nbrligne][j]/contingence2[nbrligne][nbrcolonnes];
            // //cout<<"="<<contingencetheo[i][j]<<endl;
          }
      }
    for (int i=0;i<2;i++){
      for (int j=0;j<9;j++){
        cout<<contingencetheo[i][j]<<";";
      }
      cout<<endl;
    }
    cout<<endl;
          float scorekhi2=0;
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
          //                    //cout << "Not reliable test, ";
          //                    //cout << "expected: " << e << endl;
                          test=1;
                      }
                  }
              }
            }
            for(unsigned i=0; i<nbrligne; ++i)
            {
               for(unsigned j=0; j<nbrcolonnes; ++j)
               {
                     if(contingence2[i][j] < 5)
                     test=1;
               }
            }
            //cout<<scorekhi2<<endl;
            //cout<<"test:"<<test<<endl;
              float pval;
              if(test==0)
            {           /* for(int i(0); i<(nbrligne); ++i)
                        {    for(int j(0); j<(nbrcolonnes); ++j)
                                        {
                                          scorekhi2 += (pow(((contingence2[i][j]-(contingencetheo[i][j]))),2)/contingencetheo[i][j]);
                                           //cout<<scorekhi2<<"="<<contingence2[i][j]<<"-"<<contingencetheo[i][j]<<"^2"<<"/"<<contingencetheo[i][j]<<endl;
                                        }
                                    }*/


              for(int i(0); i<(nbrligne); ++i)
              {    for(int j(0); j<(nbrcolonnes); ++j)
                  {
                    if (contingence2[i][j] != 0 ){
                      //cout<<"contingence2[i][j]:"<<contingence2[i][j]<<endl;
                      //cout<<"ontingencetheo[i][j]:"<<contingencetheo[i][j]<<endl;

                          double div=contingence2[i][j]/(double)contingencetheo[i][j];
                          //cout<<"contingence2[i][j]:"<<contingence2[i][j]<<endl;
                          //cout<<contingence2[i][j]<<"/"<<contingencetheo[i][j]<<"="<<div<<endl;
                          //cout<<contingence2[i][j]<<"*"<<log(div)<<endl;
                          scorekhi2 += contingence2[i][j] * log(div);
                          //cout<<"scorekhi2:"<<scorekhi2<<endl;
                        }

                  }
                }
                if (scorekhi2 > 0){

                boost::math::chi_squared mydist(8);
                pval = 1 - boost::math::cdf(mydist, scorekhi2);
                if(pval == 0){
                    pval = 2.0e-16;}
                cout<<"score: "<<scorekhi2<<endl;
                cout<<"p: "<<pval<<endl;}
                else{cout<<"score: "<<scorekhi2<<endl;}
            }
            else {
              scorekhi2=0;
              pval=1;
              //cout<<"score: "<<endl;
              //cout<<"p: "<<endl;
              continue;         }
              patternscore p1;
              p1.pattern1=snpNameList[l1];
              p1.pattern2=snpNameList[l2];
              p1.score=scorekhi2;
              patternscoreList.push_back(p1);

                    }
                      }

}

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
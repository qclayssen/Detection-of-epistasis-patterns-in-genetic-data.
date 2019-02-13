//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
//Statistique fonction for Epistasis detection
//Created :09/11/18
//Modified :11/02/2019

#include "../include/stats.hpp"


//=================================================
//create contingency table pattern2
//=================================================
//Goes through the genotype matrix and completes a contingency table for a pattern of size 2
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




//=================================================
//create contingency table pattern3
//=================================================
//Goes through the genotype matrix and completes a contingency table for a pattern of size 3

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



//=================================================
//g test 2SNP
//=================================================
//Calculates gscore and p-value thanks to contingency table and boost
//for a pattern with 2 snp
score_pval g_test_2SNP(contingence2SNP contingence2){
  int const nbrcolonnes(9);
  int const nbrligne(2);
  score_pval result_gtest;

  for(int j=0; j<(nbrcolonnes); ++j)
  {    for(int i=0; i<(nbrligne); ++i)
      {
         contingence2[nbrligne][j] += contingence2[i][j];

      }

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

  float contingencetheo[nbrligne][nbrcolonnes];
  for (int i=0;i<(nbrligne);i++){ //Set contingency matrix to 0
    for (int j=0;j<(nbrcolonnes);j++){
      contingencetheo[i][j]=0;
    }
  }
  for(int i=0; i<(nbrligne); ++i){
    for(int j=0; j<(nbrcolonnes); ++j){
      contingencetheo[i][j] = (float) contingence2[i][nbrcolonnes]* (float)contingence2[nbrligne][j]/(float)contingence2[nbrligne][nbrcolonnes];
    }
  }

float scorekhi2=0;
int test;
test=0;

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

  int df;
  test=0;
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
  }
  else{
    scorekhi2=0;
    pval=1;
  }
  result_gtest.score=scorekhi2;
  result_gtest.pval=pval;
  return(result_gtest);
}


//=================================================
//g test 3SNP
//=================================================
//Calculates gscore and p-value thanks to contingency table and boost
//for a pattern with 3 snp
score_pval g_test_3SNP(contingence3SNP contingence2){
  int const nbrcolonnes(27);
  int const nbrligne(2);
  score_pval result_gtest;

  for(int j=0; j<(nbrcolonnes); ++j)
  {    for(int i=0; i<(nbrligne); ++i)
      {
         contingence2[nbrligne][j] += contingence2[i][j];
      }
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

  float contingencetheo[nbrligne][nbrcolonnes];
  for (int i=0;i<(nbrligne);i++){ //Set contingency matrix to 0
    for (int j=0;j<(nbrcolonnes);j++){
      contingencetheo[i][j]=0;
    }
  }
  for(int i=0; i<(nbrligne); ++i){
    for(int j=0; j<(nbrcolonnes); ++j){
      contingencetheo[i][j] = (float) contingence2[i][nbrcolonnes]* (float)contingence2[nbrligne][j]/(float)contingence2[nbrligne][nbrcolonnes];
    }
  }

float scorekhi2=0;

int test;
test=0;

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

  int df;
  test=0;
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
  }
  else{
    scorekhi2=0;
    pval=1;
  }
  result_gtest.score=scorekhi2;
  result_gtest.pval=pval;
  return(result_gtest);
}


//=================================================
//add_gtest_pval
//=================================================
//Adds pval in the structure of solution
float add_gtest_pval (patternscore pattern,blas_matrix genos,blas_matrix phenos_m){
  score_pval bi;
  float score;
  if(pattern.snp3==-1){//for a pattern with 2 snp
    contingence2SNP contingence2;
    contingence2SNP* adr_contingence2 = &contingence2;
    create_contingency_table_pattern2(pattern.snp1,pattern.snp2,adr_contingence2,genos,phenos_m);
    bi=g_test_2SNP(contingence2);
    score=bi.pval;
  }
  else{//for a pattern with 3 snp
    contingence3SNP contingence3;
    contingence3SNP* adr_contingence3 = &contingence3;
    create_contingency_table_pattern3(pattern.snp1,pattern.snp2,pattern.snp3,adr_contingence3,genos,phenos_m);
    bi=g_test_3SNP(contingence3);
    score=bi.pval;
  }
  return(score);
}

//=================================================
//add_gtest_score
//=================================================
//Adds pval in the structure of solution
float add_gtest_score (patternscore pattern,blas_matrix genos,blas_matrix phenos_m){
  score_pval bi;
  float score;
  if(pattern.snp3==-1){
    contingence2SNP contingence2;
    contingence2SNP* adr_contingence2 = &contingence2;
    create_contingency_table_pattern2(pattern.snp1,pattern.snp2,adr_contingence2,genos,phenos_m);
    bi=g_test_2SNP(contingence2);
    score=bi.score;
  }
  else{
    contingence3SNP contingence3;
    contingence3SNP* adr_contingence3 = &contingence3;
    create_contingency_table_pattern3(pattern.snp1,pattern.snp2,pattern.snp3,adr_contingence3,genos,phenos_m);
    bi=g_test_3SNP(contingence3);
    score=bi.score;
  }
  return(score);
}


//=================================================
//add_gtest_results
//=================================================
//Returns both pval and score of a given pattern
score_pval add_gtest_results (patternscore pattern,blas_matrix genos,blas_matrix phenos_m){
  score_pval bi;
  if(pattern.snp3==-1){
    contingence2SNP contingence2;
    contingence2SNP* adr_contingence2 = &contingence2;
    create_contingency_table_pattern2(pattern.snp1,pattern.snp2,adr_contingence2,genos,phenos_m);
    bi=g_test_2SNP(contingence2);
  }
  else{
    contingence3SNP contingence3;
    contingence3SNP* adr_contingence3 = &contingence3;
    create_contingency_table_pattern3(pattern.snp1,pattern.snp2,pattern.snp3,adr_contingence3,genos,phenos_m);
    bi=g_test_3SNP(contingence3);
  }
  return(bi);
}

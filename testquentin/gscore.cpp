#include "gscore.hpp"

void MyClass::foo()
{    for (l1=0;l1<int(genos.size2())-1;l1++){ //First SNP of the pattern
      for (l2=l1+1;l2<int(genos.size2());l2++){ //Second SNP of the pattern
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

        float scorekhi2=0;/*
        for(int i(0); i<(nbrligne); ++i)
        {    for(int j(0); j<(nbrcolonnes); ++j)
            {
               scorekhi2 += (pow(((contingence2[i][j]-(contingencetheo[i][j]))),2)/contingencetheo[i][j]);
               //cout<<scorekhi2<<"="<<contingence2[i][j]<<"-"<<contingencetheo[i][j]<<"^2"<<"/"<<contingencetheo[i][j]<<endl;
            }
        }*/

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
      if(test==1){
        cout<<"The test isn't reliable"<<endl;
      }


        int df;
        df=(nbrligne-1)*(nbrcolonnes-1);
        float pval;
        if(test==0){
          for(int i(0); i<(nbrligne); ++i){
            for(int j(0); j<(nbrcolonnes); ++j){
              if (contingence2[i][j] != 0 ){
                double div=(double)contingence2[i][j]/contingencetheo[i][j];
                cout<<log(div)<<endl;
                scorekhi2 += contingence2[i][j] * log(div);
                cout<<scorekhi2<<endl;
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
          cout<<"score: "<<endl;
          cout<<"p: "<<endl;
          continue;
        }

        patternscore p1;
        p1.pattern1=snpNameList[l1];
        p1.pattern2=snpNameList[l2];
        p1.score=scorekhi2;
        patternscoreList.push_back(p1);
      }
    }
}

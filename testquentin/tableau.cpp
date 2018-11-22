/*

#Autor:Antoine Laiwnee (Master 2 Bioinformatics, University of Nantes)
#Autor:Quentin Clayssen (Master 2 Bioinformatics, University of Nantes)
#methode: path-relinking for episatasie detection
or //script for ...
#University of Nantes
#undersupervison: Christine Sinoquet
#dates

F mesure temps d'execution
*/
#include <math.h>
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



float khi2(int contingence[2][9] ,int nbrligne, int nbrcolonnes)
{

for(int j(0); j<(nbrcolonnes); ++j)
{    for(int i(0); i<(nbrligne); ++i)
    {

       cout<<contingence[nbrligne][j]<<"+"<<contingence[i][j]<<"=";
       contingence[nbrligne][j] += contingence[i][j];
       cout<<contingence[nbrligne][j]<<" ";

    }cout<<";"<<endl;
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
       cout<<contingence[i][nbrcolonnes]<<"*"<<contingence[nbrligne][j]<<"/"<<contingence[nbrligne][nbrcolonnes];
       contingencetheo[i][j] = contingence[i][nbrcolonnes]*contingence[nbrligne][j]/contingence[nbrligne][nbrcolonnes];
       cout<<"="<<contingencetheo[i][j]<<endl;
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
       scorekhi2 += (pow(((contingence[i][j]-(contingencetheo[i][j]))),2)/contingencetheo[i][j]);
       cout<<scorekhi2<<"="<<contingence[i][j]<<"-"<<contingencetheo[i][j]<<"^2"<<"/"<<contingencetheo[i][j]<<endl;
    }
}

float scorekhi2=0;
for(int i(0); i<(nbrligne); ++i)
{    for(int j(0); j<(nbrcolonnes); ++j)
    {
        if(c(i,j) != 0)
        {
            double div = (double) contingence[i][j] / contingencetheo[i][j]);
            scorekhi2 += contingence[i][j] * log(div);
        }
    }
}
scorekhi2  *= 2;


float p;
cout<<"score: "<<scorekhi2<<endl;
boost::math::chi_squared mydist(8);
pval = 1 - boost::math::cdf(mydist, scorekhi2)
p=boost::math::cdf(mydist,scorekhi2);
cout<<"p: "<<p<<endl;
return(scorekhi2);
}

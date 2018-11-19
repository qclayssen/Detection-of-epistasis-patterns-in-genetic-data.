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

#include <iostream>
using namespace std;



int const taillex(4);
int const tailley(3);
int matrice[taillex][tailley];

double khi(double matrice[][] ,int taillex, int tailley)
{/*
  int matrice[taillex][tailley];
     for(int i(0); i<taillex-1; ++i)
     {
        matrice[tailley][0] += matrice[i][0];
     }
  int matrice[tailley][1];
    for(int i(0); i<taillex-1; ++i)
    {
       matrice[tailley][1] += matrice[i][1];
    }
  int matrice[0][4];
    for(int i(0); i<tailley-1; ++i)
    {
       matrice[0][4] += matrice[0][i];
    }
  int matrice[1][4];
    for(int i(0); i<tailley-1; ++i)
    {
       matrice[1][4] += matrice[1][i];
    }
  int matrice[2][4];
    for(int i(0); i<tailley-1; ++i)
    {
       matrice[2][4] += matrice[2][i];
    }
  int matrice[3][4];
    for(int i(0); i<tailley-1; ++i)
    {
       matrice[3][4] += matrice[i][3];
    }*/


    for(int i(0); i<tailley-1; ++i)
    {    for(int j(0); j<taillex-1; ++i)
        {
           matrice[taillex][i] += matrice[i][j];
        }
    }


    for(int i(0); i<tailley-1; ++i)
    {    for(int j(0); j<taillex-1; ++i)
        {
           matrice[j][tailley] += matrice[i][j];
        }
    }



  int matricetheo[tailley-1][taillex-1];
    for(int i(0); i<tailley-1; ++i)
    {    for(int j(0); j<taillex-1; ++i)
        {
           matricetheo[i][j] = matrice[i][taillex]*matrice[tailley][j]/matrice[tailley][taillex];
        }
    }


    for(int i(0); i<tailley-1; ++i)
    {    for(int j(0); j<taillex-1; ++i)
        {
           matrice[taillex][i] += matrice[i][x];
        }
    }

  int scorekhi2=0;
    for(int i(0); i<tailley-1; ++i)
    {    for(int j(0); j<taillex-1; ++i)
        {
           scorekhi2 += (matrice[i][j]-pow(matricetheo[i][j],2)/matrice[tailley][taillex];
        }
    }

    cout<<scorekhi2
}

/*   double sommeattein(0);
   for(int i(0); i<taillematricex; ++i)
   {
      sommemalade += matrice[i][1];   //On additionne toutes les valeurs
   }
   return sommemalade;

   double sommenonattein(0);
   for(int i(0); i<taillematricex; ++i)
   {
      sommenonattein += matrice[i][2];   //On additionne toutes les valeurs
   }
   return sommenonattein;

   double somme0(0);
   for(int i(0); i<taillematricex; ++i)
   {
      somme0 += matrice[1][i];   //On additionne toutes les valeurs
   }
   return sommenonattein; */

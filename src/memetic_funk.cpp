//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
//Fonction of memetic algorithm for Epistasis detection
//Created :09/11/18
//Modified :11/02/2019


#include "../include/memetic_funk.hpp"



//=================================================
// Initialisation of the population
//=================================================


// return a solution vector of size n randomly drawn from the pattern list
vector<patternscore> initialize_population(int n,vector<patternscore>* adr_patternscoreList){
  vector<patternscore> pop;
  srand(time(0));//Generation of a new random seed
  random_shuffle((*adr_patternscoreList).begin(), (*adr_patternscoreList).end());//Randomly shuffles all the solutions
  for (unsigned int i=0;i<n;i++){
      pop.push_back((*adr_patternscoreList)[i]);
  }
  return(pop);
}



//=================================================
// select_pairs_of_individuals_to_be_crossed
//=================================================


// create a couple vector of size n_p randomly drawn from the population, parents are following each other in the randomized population list
vector<parents_pairs> select_pairs_of_individuals_to_be_crossed(vector<patternscore> n_pairs_selected_parents,int n_p){
  vector<parents_pairs>pairs_of_parents;
  parents_pairs couple;
  srand(time(0));
  random_shuffle(n_pairs_selected_parents.begin(), n_pairs_selected_parents.end());
  for (int i=0;i<n_p;i+=2){//the loop advances by couple
    if(i+2>n_pairs_selected_parents.size() || i+1>n_pairs_selected_parents.size()){
      break;
    }
    couple.parent1=n_pairs_selected_parents[i];
    couple.parent2=n_pairs_selected_parents[i+1];
    pairs_of_parents.push_back(couple);
  }return(pairs_of_parents);
}



//=================================================
// create_two_children_for_each_selected_pair_of_parents
//=================================================


// for each couple we creat 2 solution (son and daughter)

vector<patternscore> create_two_children_for_each_selected_pair_of_parents(vector<parents_pairs> pairs_of_parents){
   patternscore son;
   patternscore daugther;

   vector<patternscore> children_parents;
   for (int i=0;i<pairs_of_parents.size();++i){

     if(pairs_of_parents[i].parent1.snp1==pairs_of_parents[i].parent2.snp1) //the son inherits from the father and an snp from the mother
       {son=pairs_of_parents[i].parent1;
        son.snp2=pairs_of_parents[i].parent2.snp2;
        children_parents.push_back(son);}
     son=pairs_of_parents[i].parent1;
     son.snp2=pairs_of_parents[i].parent2.snp1;
     children_parents.push_back(son);

     if(pairs_of_parents[i].parent1.snp2==pairs_of_parents[i].parent2.snp2)//the daughter inherits from the mother and an snp from the father
      {daugther=pairs_of_parents[i].parent2;
           daugther.snp2=pairs_of_parents[i].parent1.snp1;
           children_parents.push_back(daugther);}
     daugther=pairs_of_parents[i].parent2;
     daugther.snp2=pairs_of_parents[i].parent1.snp2;
     children_parents.push_back(daugther);

   } return(children_parents);


 }



//=================================================
// perform_one_mutation_per_child: mutation of child base on random
//=================================================

// Change a snp of child randomly
void perform_one_mutation_per_child(vector<patternscore>* adr_children_parents,int prob_mutation,vector<patternscore> patternscoreList){
 for (int i=0;i<(*adr_children_parents).size();++i){
   int mutation = rand() % 100 + 0;
   if (mutation >= prob_mutation){
     break;
   }else{
     int snp = rand() % patternscoreList.size();// Select snp donor
     int parentpattern;
     int mutpattern;
     if (patternscoreList[snp].snp3 != -1)//snp select  from donor
       {parentpattern = rand() % 3;}
     else
       {parentpattern = rand() % 2;}

     if ((*adr_children_parents)[i].snp3 != -1)//snp that's going to be changed
       {mutpattern = rand() % 3 ;}
     else
       {mutpattern = rand() % 2;}

       switch (parentpattern) {
         case 0:
          switch (mutpattern) {
            case 0:
            (*adr_children_parents)[i].snp1=patternscoreList[snp].snp1;
            break;
            case 1:
            (*adr_children_parents)[i].snp1=patternscoreList[snp].snp2;
            break;
            case 2:
            (*adr_children_parents)[i].snp1=patternscoreList[snp].snp3;
            break;}
         case 1:
           switch (mutpattern) {
             case 0:
             (*adr_children_parents)[i].snp2=patternscoreList[snp].snp1;
             break;
             case 1:
             (*adr_children_parents)[i].snp2=patternscoreList[snp].snp2;
             break;
             case 2:
             (*adr_children_parents)[i].snp2=patternscoreList[snp].snp3;
             break;}
         case 2:
           switch (mutpattern) {
            case 0:
            (*adr_children_parents)[i].snp3=patternscoreList[snp].snp1;
            break;
            case 1:
            (*adr_children_parents)[i].snp3=patternscoreList[snp].snp2;
            break;
            case 2:
            (*adr_children_parents)[i].snp3=patternscoreList[snp].snp3;
            break;}
         }
     }
   }
 }



//=================================================
//Sort of solutions
//=================================================


vector<patternscore> identify_best_solutions(vector<patternscore> pop, int k, int n){
  vector<patternscore> best_solutions;
  if (pop[1].pval!=0){
    std::sort(pop.begin(), pop.end(), compareByScore);
  }
  else{
    std::sort(pop.begin(), pop.end(), compareByPval);
  }
  int i=0;
  while(i<k){
        best_solutions.push_back(pop[i]);
        i=i+1;}
  return(best_solutions);
  }


//=================================================
//update population
//=================================================

//replace parent by there children if they have netter score/pvalue

void update_population(vector<patternscore> children_parents, vector<patternscore>* adr_pop,int n){
      for (int j=0;j<(*adr_pop).size();j++){
        for (int i=0;i<children_parents.size();i++){
          if((*adr_pop)[j].idparent==children_parents[i].idparent){
            if ((*adr_pop)[j].pval!=0){
              if (children_parents[i].pval>(*adr_pop)[j].pval){
                (*adr_pop)[j]=children_parents[i];
              }
            }
            else{
              if (children_parents[i].score>(*adr_pop)[j].score)
              {(*adr_pop)[j]=children_parents[i];}
            }
          }
        }
      }
    }

//Autors : Quentin Clayssen, Antoine Laine (Master2 Bioinformatics, University of Nantes)
//Parsing of CSV Files for Epistasis detection
//Created :09/11/18
//Modified :11/02/2019

#include "../include/path_relinking_func.hpp"

//=================================================
/*
Initialisation of the elite solutions : they are selected at random among all possible solutions.
Returns a vector of k solutions, with k being set in the Parameters file.
*/
//=================================================
vector<patternscore> initialize_elite_solutions(unsigned int k,vector<patternscore> patternscoreList){
  vector<patternscore> elite_sols;
  vector<patternscore> nodes(patternscoreList.begin(),patternscoreList.end());
  srand(time(0)); //Generation of a new random seed
  random_shuffle(nodes.begin(), nodes.end()); //Randomly shuffles all the solutions
  for (vector<patternscore>::iterator it=nodes.begin(); it!=nodes.end(); ++it){ //Get the k first solutions of the shuffled list
    if(elite_sols.size()<k){
      elite_sols.push_back(*it);
    }
    else{
      break;
    }
  }
  return(elite_sols);
}


//=================================================
/*
Select the two solutions sA and sB : they are selected at random among the elite solutions.
Returns a vector of 2 solutions, the first one considered as sA, the second one considered as sB.
*/
//=================================================
vector<patternscore> select_two_solutions_at_random(vector<patternscore> elite_sols){
  vector<patternscore> sA_sB;
  vector<patternscore> nodes(elite_sols.begin(),elite_sols.end());
  srand(time(0)); //Generation of a new random seed
  random_shuffle(nodes.begin(), nodes.end()); //Randomly shuffles the elite solutions
  if(nodes[0].pval > nodes[1].pval){ //The solutions with the smallest p-value is pushed first in the returned vector (and thus, it becomes sA)
    sA_sB.push_back(nodes[0]);
    sA_sB.push_back(nodes[1]);
  }
  else{
    sA_sB.push_back(nodes[1]);
    sA_sB.push_back(nodes[0]);
  }
  return(sA_sB);
}

//=================================================
/*
Select in the neighbourhood of sA the solution that is the closest to sB.
Returns a unique solution.
*/
//=================================================
patternscore select_closest_neighbor_to_guiding_solution(patternscore s,patternscore sB, vector<patternscore> patternscoreList,int s_n){
  patternscore s_closest_neighbour;
  vector<patternscore> s_neighbours=neighbours(s,patternscoreList,s_n); //Create a vector containing s_n neighbours of sA
  int min_delta=999;
  for (unsigned int i=0;i<s_neighbours.size();i++){ //In the vector previously created, searches the solution with the minimal delta with sB and returns it.
    int delta=calculate_delta(s_neighbours[i],sB);
    if (delta<min_delta){
      min_delta=delta;
      s_closest_neighbour=s_neighbours[i];
    }
  }
  return(s_closest_neighbour);
}

//=================================================
/*
Check if the p-value/score of the input solution is better than any of the values in the elite solutions
Returns 1 if the solution is better than an elite solution, returns 0 if it isn't.
*/
//=================================================
int promizing_score(patternscore s_closest_neighbour,vector<patternscore> elite_sols){
  int min_elite=0;
  double min_score_pval=0;
  double min_score=999999;
  int score_or_pval=0;
  for (unsigned int i=0;i<elite_sols.size();i++){ //Checks if p-value = 0 (due to datas), if so, it will compare the score instead.
    if (elite_sols[i].pval==0){
      score_or_pval=1;
    }
  }
  if (score_or_pval==0){ //Comparing P-value
    for (unsigned int i=0;i<elite_sols.size();i++){ //Compare the solution with every elite solution, and returns 1 if it's find better than at least one.
      if (elite_sols[i].pval>min_score_pval){
        min_score_pval=elite_sols[i].pval;
        min_elite=i;
      }
    }
    if (elite_sols[min_elite].pval>s_closest_neighbour.pval){
      return(1);
    }
    else{
      return(0);
    }
  }
  else{ //Comparing scores
    for (unsigned int i=0;i<elite_sols.size();i++){ //Compare the solution with every elite solution, and returns 1 if it's find better than at least one.
      if (elite_sols[i].score<min_score){
        min_score=elite_sols[i].score;
        min_elite=i;
      }
    }
    if (elite_sols[min_elite].score<s_closest_neighbour.score){
      return(1);
    }
    else{
      return(0);
    }
  }
}

//=================================================
/*
Replace the value with the lowest score among the elite solution by an optimized solution.
The vector of elite solutions is given by adress, so there is no return.
*/
//=================================================
void update(patternscore s_opt, vector<patternscore>* adr_elite_sols){
  double min_score_pval=0;
  double min_score=999999;
  int min_elite=0;
  int score_or_pval=0;
  for (unsigned int i=0;i<(*adr_elite_sols).size();i++){ //Checks if p-value = 0 (due to datas), if so, it will compare the score instead.
    if ((*adr_elite_sols)[i].pval==0){
      score_or_pval=1;
    }
  }
  if (score_or_pval==0){ //Compare p-value
    for(unsigned int i=0;i<(*adr_elite_sols).size();i++){ //Compare all the elite sols to find the one with the lowest score
      if ((*adr_elite_sols)[i].snp1==s_opt.snp1 && (*adr_elite_sols)[i].snp2==s_opt.snp2){
        return;
      }
      if ((*adr_elite_sols)[i].pval>min_score_pval){
        min_score_pval=(*adr_elite_sols)[i].pval;
        min_elite=i;
      }
    }
    (*adr_elite_sols)[min_elite]=s_opt; //Replace the lowest solution by the optimized solution.
  }
  else{ //Compare scores
    for(unsigned int i=0;i<(*adr_elite_sols).size();i++){ //Compare all the elite sols to find the one with the lowest score
      if ((*adr_elite_sols)[i].snp1==s_opt.snp1 && (*adr_elite_sols)[i].snp2==s_opt.snp2){
        return;
      }
      if ((*adr_elite_sols)[i].score<min_score){
        min_score=(*adr_elite_sols)[i].score;
        min_elite=i;
      }
    }
    (*adr_elite_sols)[min_elite]=s_opt; //Replace the lowest solution by the optimized solution.
  }
}

//=================================================
/*
Out File Manager, formating the results and writing them in a file located in the outputs directory.
On this file, the elite solutions with there scores and p-value, plus the time that was necessary for this execution of the method.
*/
//=================================================
void outfilePR(string genos_file ,vector<string> snpNameList,vector<patternscore> best_solutions,int s_n, float duree){
  //The name of the file is generated with the name of the genotype file that was in input of the method.
  string file_basename = basename((char*)genos_file.c_str());
  string result_filename = "outputs/s_n"+to_string(s_n)+"_"+file_basename;
  std::ofstream _results_handler;
  _results_handler.open(result_filename.c_str(), ios::trunc);

  if(!_results_handler)
  {
      std::cerr << "Error while opening output.txt (by writing access) !\n";
      exit(-1);
  }

  if (_results_handler.is_open()){ //If the opening of the file is a success, the results are written in a set format.
      _results_handler<<"Pattern"<<"\t"<<"p-value"<<"\t"<<"score"<<endl;
      for (vector<patternscore>::iterator it=best_solutions.begin();it!=best_solutions.end();it++){
        if ((*it).snp3==-1){
          _results_handler<<"<"<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<">"<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl<<"durée: "<<duree<<" seconde"<<endl;
        }
        else{
          _results_handler<<"<"<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<","<<snpNameList[(*it).snp3]<<">"<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl;
        }
      }_results_handler<<endl<<"durée: "<<duree<<" seconds"<<endl;
    }
  else
  {
    cout<<"fail"<<endl;
  }
}

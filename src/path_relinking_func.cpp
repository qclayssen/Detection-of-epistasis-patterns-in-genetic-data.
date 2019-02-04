#include "../include/path_relinking_func.hpp"


vector<patternscore> initialize_elite_solutions(unsigned int k,vector<patternscore> patternscoreList){
  vector<patternscore> elite_sols;
  vector<patternscore> nodes(patternscoreList.begin(),patternscoreList.end());
  srand(time(0));
  random_shuffle(nodes.begin(), nodes.end());
  for (vector<patternscore>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
    if(elite_sols.size()<k){
      elite_sols.push_back(*it);
    }
    else{
      break;
    }
  }
  return(elite_sols);
}

vector<patternscore> select_two_solutions_at_random(vector<patternscore> elite_sols){
  vector<patternscore> sA_sB;
  vector<patternscore> nodes(elite_sols.begin(),elite_sols.end());
  srand(time(0));
  random_shuffle(nodes.begin(), nodes.end());
  if(nodes[0].pval > nodes[1].pval){
    sA_sB.push_back(nodes[0]);
    sA_sB.push_back(nodes[1]);
  }
  else{
    sA_sB.push_back(nodes[1]);
    sA_sB.push_back(nodes[0]);
  }
  return(sA_sB);
}

patternscore select_closest_neighbor_to_guiding_solution(patternscore s,patternscore sB, vector<patternscore> patternscoreList, int s_n){
  patternscore s_closest_neighbour;
  vector<patternscore> s_neighbours=neighbours(s,patternscoreList,s_n);
  int min_delta=999;
  for (unsigned int i=0;i<s_neighbours.size();i++){
    int delta=calculate_delta(s_neighbours[i],sB);
    if (delta<min_delta){
      min_delta=delta;
      s_closest_neighbour=s_neighbours[i];
    }
  }
  return(s_closest_neighbour);
}

int promizing_score(patternscore s_closest_neighbour,vector<patternscore> elite_sols){
  int min_elite;
  double min_score_pval=0;
  double min_score=999999;
  int score_or_pval=0;
  for (unsigned int i=0;i<elite_sols.size();i++){
    if (elite_sols[i].pval==0){
      score_or_pval=1;
    }
  }
  if (score_or_pval=0){
    for (unsigned int i=0;i<elite_sols.size();i++){
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
  else{
    for (unsigned int i=0;i<elite_sols.size();i++){
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


void update(patternscore s_opt, vector<patternscore>* adr_elite_sols){
  double min_score_pval=0;
  double min_score=999999;
  int min_elite=0;
  int score_or_pval=0;
  for (unsigned int i=0;i<(*adr_elite_sols).size();i++){
    if ((*adr_elite_sols)[i].pval==0){
      score_or_pval=1;
    }
  }
  if (score_or_pval=0){
    for(unsigned int i=0;i<(*adr_elite_sols).size();i++){
      if ((*adr_elite_sols)[i].snp1==s_opt.snp1 && (*adr_elite_sols)[i].snp2==s_opt.snp2){
        return;
      }
      if ((*adr_elite_sols)[i].pval>min_score_pval){
        min_score_pval=(*adr_elite_sols)[i].pval;
        min_elite=i;
      }
    }
    (*adr_elite_sols)[min_elite]=s_opt;
  }
  else{
    for(unsigned int i=0;i<(*adr_elite_sols).size();i++){
      if ((*adr_elite_sols)[i].snp1==s_opt.snp1 && (*adr_elite_sols)[i].snp2==s_opt.snp2){
        return;
      }
      if ((*adr_elite_sols)[i].score<min_score){
        min_score=(*adr_elite_sols)[i].score;
        min_elite=i;
      }
    }
    (*adr_elite_sols)[min_elite]=s_opt;
  }
}

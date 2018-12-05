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
  if(nodes[0].score < nodes[1].score){
    sA_sB.push_back(nodes[0]);
    sA_sB.push_back(nodes[1]);
  }
  else{
    sA_sB.push_back(nodes[1]);
    sA_sB.push_back(nodes[0]);
  }
  return(sA_sB);
}


int calculate_delta(patternscore s, patternscore sB){
  int diff=0;
  if (s.pattern1!=sB.pattern1 && s.pattern1!=sB.pattern2 && s.pattern1!=sB.pattern3){
    diff=diff+1;
  }
  if (s.pattern2!=sB.pattern1 && s.pattern2!=sB.pattern2 && s.pattern2!=sB.pattern3){
    diff=diff+1;
  }
  if (s.pattern3!=sB.pattern1 && s.pattern3!=sB.pattern2 && s.pattern3!=sB.pattern3){
    diff=diff+1;
  }
  return(diff);
}

vector<patternscore> neighbours(patternscore s,vector<patternscore> patternscoreList){
  vector<patternscore> s_neighbours;
  for (unsigned int i=0;i<patternscoreList.size();i++){
    int delta=calculate_delta(s,patternscoreList[i]);
    if (delta==1){
      s_neighbours.push_back(patternscoreList[i]);
    }
  }
  return(s_neighbours);
}

patternscore select_closest_neighbor_to_guiding_solution(patternscore s,patternscore sB, vector<patternscore> patternscoreList){
  patternscore s_closest_neighbour;
  vector<patternscore> s_neighbours=neighbours(s,patternscoreList);
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
  double min_score=99999999;
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

patternscore hill_climbing_lc(patternscore s_closest_neighbour, vector<patternscore> patternscoreList,blas_matrix genos,blas_matrix phenos_m){
  vector<patternscore> s_neighbours = neighbours(s_closest_neighbour,patternscoreList);
  patternscore actual_s = s_closest_neighbour;
  for (unsigned int i=0;i<s_neighbours.size();i++){
    s_neighbours[i].score=add_gtest_score(s_neighbours[i],genos,phenos_m);
    if (s_neighbours[i].score>actual_s.score){
      actual_s=s_neighbours[i];
    }
  }
  return(actual_s);
}

void update(patternscore s_opt, vector<patternscore>* adr_elite_sols){
  double min_score=99999999;
  int min_elite=0;
  for(unsigned int i=0;i<(*adr_elite_sols).size();i++){
    if ((*adr_elite_sols)[i].pattern1==s_opt.pattern1 && (*adr_elite_sols)[i].pattern2==s_opt.pattern2){
      return;
    }
    if ((*adr_elite_sols)[i].score<min_score){
      min_score=(*adr_elite_sols)[i].score;
      min_elite=i;
    }
  }
  (*adr_elite_sols)[min_elite]=s_opt;
}

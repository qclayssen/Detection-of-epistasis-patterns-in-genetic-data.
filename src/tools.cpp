#include "../include/tools.hpp"


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

void cout_list(vector<patternscore> list_to_cout,vector<string> snpNameList){
  for (vector<patternscore>::iterator it=list_to_cout.begin();it!=list_to_cout.end();it++){
    if ((*it).snp3==-1){
      cout<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<endl;
      cout<<(*it).score<<endl;
    }
    else{
      cout<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<","<<snpNameList[(*it).snp3]<<endl;
      cout<<(*it).score<<endl;
    }
  }
}

int calculate_delta(patternscore s, patternscore sB){
  int diff=0;
  if (s.snp1!=sB.snp1 && s.snp1!=sB.snp2 && s.snp1!=sB.snp3){
    diff=diff+1;
  }
  if (s.snp2!=sB.snp1 && s.snp2!=sB.snp2 && s.snp2!=sB.snp3){
    diff=diff+1;
  }
  if (s.snp3!=sB.snp1 && s.snp3!=sB.snp2 && s.snp3!=sB.snp3){
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

patternscore hill_climbing_lc(patternscore s_closest_neighbour, vector<patternscore> patternscoreList,blas_matrix genos,blas_matrix phenos_m){
  vector<patternscore> s_neighbours = neighbours(s_closest_neighbour,patternscoreList);
  patternscore actual_s = s_closest_neighbour;
  for (unsigned int i=0;i<s_neighbours.size();i++){
    s_neighbours[i].score=add_gtest_score(s_neighbours[i],genos,phenos_m);
    if (s_neighbours[i].score<actual_s.score){
      actual_s=s_neighbours[i];
    }
  }
  return(actual_s);
}


void outfile(string filename, vector<patternscore> list_to_cout,vector<string> snpNameList){  ofstream filename;
  filename.open ("example.txt");
  filename << cout_list(best_solutions,snpNameList);
  filename.close();
  return 0;}

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
      cout<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl;
    }
    else{
      cout<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<","<<snpNameList[(*it).snp3]<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl;
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
  srand(time(0));
  random_shuffle(patternscoreList.begin(), patternscoreList.end());
  int i=0;
  int j=0;
  while (j<2){
    int delta=calculate_delta(s,patternscoreList[i]);
    if (delta==1){
      s_neighbours.push_back(patternscoreList[i]);
      j++;
    }
    cout<<i<<endl;
    i++;
  }
  return(s_neighbours);
}

patternscore hill_climbing_lc(patternscore s_closest_neighbour, vector<patternscore> patternscoreList,blas_matrix genos,blas_matrix phenos_m,int s_n){
  /*vector<patternscore> s_neighbours2 = neighbours(s_closest_neighbour,patternscoreList);
  srand(time(0));
  random_shuffle(s_neighbours2.begin(), s_neighbours2.end());
  vector<patternscore> s_neighbours;
  for (unsigned int i=0;i<s_n;i++){
    s_neighbours.push_back(s_neighbours2[i]);
  }*/
  auto taa = Clock::now();
  vector<patternscore> s_neighbours = neighbours(s_closest_neighbour,patternscoreList);
  auto tbb = Clock::now();
  patternscore actual_s = s_closest_neighbour;
  score_pval biScore;
  for (unsigned int i=0;i<s_neighbours.size();i++){
    biScore=add_gtest_results(s_neighbours[i],genos,phenos_m);
    s_neighbours[i].score=biScore.score;
    s_neighbours[i].pval=biScore.pval;
    if (s_neighbours[i].pval!=0){
      if (s_neighbours[i].pval<actual_s.pval){
        actual_s=s_neighbours[i];
      }
    }
    else{
      if (s_neighbours[i].score>actual_s.score){
        actual_s=s_neighbours[i];
      }
    }
  }
  auto tcc = Clock::now();
  std::cout << "before:"
            << duration_cast<duration<double>>(tbb - taa).count()
            << " seconds" << std::endl;
  std::cout << "after:"
            << duration_cast<duration<double>>(tcc - tbb).count()
            << " seconds" << std::endl;
  return(actual_s);
}

void outfile(string genos_file ,vector<string> snpNameList,vector<patternscore> best_solutions){
  string file_basename = basename((char*)genos_file.c_str());
  string result_filename = "outputs/RESULT_" + file_basename;
  std::ofstream _results_handler;
  _results_handler.open(result_filename.c_str(), ios::trunc);

  if(!_results_handler)
  {
      std::cerr << "Error while opening output.txt (by writing access) !\n";
      exit(-1);
  }

  if (_results_handler.is_open()){
      _results_handler<<"Pattern"<<"\t"<<"p-value"<<"\t"<<"score"<<endl;
      for (vector<patternscore>::iterator it=best_solutions.begin();it!=best_solutions.end();it++){
        if ((*it).snp3==-1){
          _results_handler<<"<"<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<">"<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl;
        }
        else{
          _results_handler<<"<"<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<","<<snpNameList[(*it).snp3]<<">"<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl;
        }
      }
    }
  else
  {
    cout<<"fail"<<endl;
  }
}

bool compareByPval(const patternscore &a, const patternscore &b){
  return a.pval < b.pval;
}

bool compareByScore(const patternscore &a, const patternscore &b){
  return a.score > b.score;
}

vector<patternscore> sort_solutions(vector<patternscore> solutions){
  int score_or_pval=0;
  for (unsigned int i=0;i<solutions.size();i++){
    if (solutions[i].pval==0){
      score_or_pval=1;
    }
  }
  if (score_or_pval==0){
    std::sort(solutions.begin(), solutions.end(), compareByPval);
  }
  else{
    std::sort(solutions.begin(), solutions.end(), compareByScore);
  }
  return(solutions);
}

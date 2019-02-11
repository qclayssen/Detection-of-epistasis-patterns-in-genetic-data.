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

vector<patternscore> neighbours(patternscore s,vector<patternscore> patternscoreList,int s_n){
  vector<patternscore> s_neighbours;
  int countForBreak=0;
  srand(time(0));
  random_shuffle(patternscoreList.begin(), patternscoreList.end());
  for (unsigned int i=0 ; i<patternscoreList.size() ; i++){
    int delta=calculate_delta(s,patternscoreList[i]);
    if (delta==1){
        s_neighbours.push_back(patternscoreList[i]);
        countForBreak+=1;
    }
    if (countForBreak==s_n){
      break;
    }
  }
  return(s_neighbours);
}

patternscore hill_climbing_lc(patternscore s_closest_neighbour, vector<patternscore> patternscoreList,blas_matrix genos,blas_matrix phenos_m,int s_n){
  vector<patternscore> s_neighbours = neighbours(s_closest_neighbour,patternscoreList,s_n);
/*  srand(time(0));
  random_shuffle(s_neighbours2.begin(), s_neighbours2.end());
  vector<patternscore> s_neighbours;
  for (unsigned int i=0;i<s_n;i++){
    s_neighbours.push_back(s_neighbours2[i]);
  }*/
  patternscore actual_s = s_closest_neighbour;
  /*auto h2 = Clock::now();
  std::cout << "neighbours: "
            << duration_cast<duration<double>>(h2 - h1).count()
            << " seconds" << std::endl;*/
  for (unsigned int i=0;i<s_neighbours.size();i++){
    score_pval biScore=add_gtest_results(s_neighbours[i],genos,phenos_m);
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
  /*auto h3 = Clock::now();
  std::cout << "actual_s: "
            << duration_cast<duration<double>>(h3 - h2).count()
            << " seconds" << std::endl;*/

  return(actual_s);
}


void outfile(string genos_file ,vector<string> snpNameList,vector<patternscore> best_solutions,int s_n, int n, float duree ,int n_it){
  string file_basename = basename((char*)genos_file.c_str());
  string result_filename = "outputs/s_n"+to_string(s_n) +"_n"+to_string(n)+"_n_it_"+to_string(n_it)+"_"+file_basename;
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
          _results_handler<<"<"<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<">"<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl<<"durée: "<<duree<<" seconde"<<endl;
        }
        else{
          _results_handler<<"<"<<snpNameList[(*it).snp1]<<","<<snpNameList[(*it).snp2]<<","<<snpNameList[(*it).snp3]<<">"<<"\t"<<(*it).pval<<"\t"<<(*it).score<<endl;
        }
      }
      _results_handler<<endl<<"durée: "<<duree<<" seconds"<<endl;
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

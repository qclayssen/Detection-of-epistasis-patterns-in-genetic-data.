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

void cout_list(vector<patternscore> list_to_cout){
  for (vector<patternscore>::iterator it=list_to_cout.begin();it!=list_to_cout.end();it++){
    if ((*it).snp3==-1){
      cout<<(*it).pattern1<<","<<(*it).pattern2<<endl;
      cout<<(*it).score<<endl;
    }
    else{
      cout<<(*it).pattern1<<","<<(*it).pattern2<<","<<(*it).pattern3<<endl;
      cout<<(*it).score<<endl;
    }
  }
}

#include "../include/memetic_funk.hpp"


/*

void cout_list(vector<patternscore> list_to_cout){
  for (vector<patternscore>::iterator it=list_to_cout.begin();it!=list_to_cout.end();it++){
    if ((*it).pattern3==""){
      hill_climbing_lc2<<(*it).snp1<<","<<(*it).snp2<<endl;
      cout<<"idparent :"<<(*it).idparent<<" score:"<<(*it).pval<<endl;
      cout<<endl;
    }
    else{
      cout<<(*it).snp1<<","<<(*it).snp2<<(*it).pattern3<<endl;

    }
  }
}
*/
vector<patternscore> initialize_population(int n,vector<patternscore> patternscoreList){
  vector<patternscore> pop;
  vector<patternscore> nodes(patternscoreList.begin(),patternscoreList.end());
  int i=0;
  srand(time(0));
  random_shuffle(nodes.begin(), nodes.end());
  for (vector<patternscore>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
    if(pop.size()<n){
      (*it).idparent=i;
      pop.push_back(*it);
      i=i+1;
    }
    else{
      break;
    }
  }
  return(pop);
}

vector<parents_pairs> select_pairs_of_individuals_to_be_crossed(vector<patternscore> n_pairs_selected_parents){
  parents_pairs couple;
  vector<parents_pairs>pairs_of_parents;
  vector<patternscore>nodes(n_pairs_selected_parents.begin(),n_pairs_selected_parents.end());
  srand(time(0));
  random_shuffle(nodes.begin(), nodes.end());
  for (int i=0;i<n_pairs_selected_parents.size();i+=2){
    if(i+2>n_pairs_selected_parents.size() || i+1>n_pairs_selected_parents.size()){
      break;
    }
    couple.parent1=nodes[i];
    couple.parent2=nodes[i+1];
    pairs_of_parents.push_back(couple);
  }return(pairs_of_parents);
}

void perform_one_mutation_per_child(vector<patternscore>* adr_children_parents,int prob_mutation){
  for (int i=0;i<(*adr_children_parents).size();++i){
    int mutation = rand() % 100;
    //cout<<mutation<<endl;
    if (mutation > prob_mutation){
      break;
      //cout<<"sheit"<<endl;
    }else{
      //int type = rand() % 3;
      int parentpattern;
      int mutpattern;
      if ((*adr_children_parents)[i].snp3 != NULL)
        {parentpattern = rand() % 3 ;}
      else
        {parentpattern = rand() % 2;}

      if ((*adr_children_parents)[i].snp3 != NULL)
        {mutpattern = rand() % 3 ;}
      else
        {mutpattern = rand() % 2;}
          int snp = rand() % (*adr_children_parents).size();
          //cout<<(*adr_children_parents)[i].snp1<<","<<(*adr_children_parents)[snp].snp1<<endl;



      if((*adr_children_parents)[i].snp1==(*adr_children_parents)[snp].snp1 || (*adr_children_parents)[i].snp2==(*adr_children_parents)[snp].snp2 || (*adr_children_parents)[i].snp3==(*adr_children_parents)[snp].snp3 )
        {continue;}

        switch (parentpattern) {
          case 1:
           switch (mutpattern) {
             case 0:
             (*adr_children_parents)[i].snp1=(*adr_children_parents)[snp].snp1;
             break;
             case 1:
             (*adr_children_parents)[i].snp1=(*adr_children_parents)[snp].snp2;
             break;
             case 2:
             (*adr_children_parents)[i].snp1=(*adr_children_parents)[snp].snp3;
             break;}
          case 2:
            switch (mutpattern) {
              case 0:
              (*adr_children_parents)[i].snp2=(*adr_children_parents)[snp].snp1;
              break;
              case 1:
              (*adr_children_parents)[i].snp2=(*adr_children_parents)[snp].snp2;
              break;
              case 2:
              (*adr_children_parents)[i].snp2=(*adr_children_parents)[snp].snp3;
              break;}
          case 3:
            switch (mutpattern) {
             case 0:
             (*adr_children_parents)[i].snp3=(*adr_children_parents)[snp].snp1;
             break;
             case 1:
             (*adr_children_parents)[i].snp3=(*adr_children_parents)[snp].snp2;
             break;
             case 2:
             (*adr_children_parents)[i].snp3=(*adr_children_parents)[snp].snp3;
             break;}
          }

            /*
          if(parentpattern == 1 ){if(mutpattern == 0){
            (*adr_children_parents)[i].snp1=(*adr_children_parents)[snp].snp1;}
            if(mutpattern == 1){
              (*adr_children_parents)[i].snp1=(*adr_children_parents)[snp].snp2;}
            if(mutpattern == 2){
              (*adr_children_parents)[i].snp1=(*adr_children_parents)[snp].snp3;}}
          if(parentpattern == 2 ){
            if(mutpattern ==0){
            (*adr_children_parents)[i].snp2=(*adr_children_parents)[snp].snp1;}
            if(mutpattern == 1){
              (*adr_children_parents)[i].snp2=(*adr_children_parents)[snp].snp2;}
            if(mutpattern == 2){
              (*adr_children_parents)[i].snp2=(*adr_children_parents)[snp].snp3;}}
          if(parentpattern == 3 ){
            if(mutpattern ==0){
            (*adr_children_parents)[i].snp3=(*adr_children_parents)[snp].snp1;}
            if(mutpattern == 1){
              (*adr_children_parents)[i].snp3=(*adr_children_parents)[snp].snp2;}
            if(mutpattern == 2){
              (*adr_children_parents)[i].snp3=(*adr_children_parents)[snp].snp3;}}*/
      }
    }
  }






vector<patternscore> create_two_children_for_each_selected_pair_of_parents(vector<parents_pairs> pairs_of_parents){
   patternscore son;
   patternscore daugther;
   vector<patternscore> children_parents;
   for (int i=0;i<pairs_of_parents.size();++i){
     if(pairs_of_parents[i].parent1.snp1==pairs_of_parents[i].parent2.snp1){continue;}
     son=pairs_of_parents[i].parent1;
     son.snp2=pairs_of_parents[i].parent2.snp1;
     children_parents.push_back(son);
     if(pairs_of_parents[i].parent1.snp2==pairs_of_parents[i].parent2.snp2){continue;}
     daugther=pairs_of_parents[i].parent2;
     daugther.snp2=pairs_of_parents[i].parent2.snp2;
     children_parents.push_back(daugther);
   } return(children_parents);


 }

/*
float add_gtest_pval (patternscore pattern,blas_matrix genos,blas_matrix phenos_m){
  float score;
  if(pattern.snp3==NULL){
    contingence2SNP contingence2;
    contingence2SNP* adr_contingence2 = &contingence2;
    cout<<pattern.pattern1<<", "<<pattern.pattern2<<endl;
    cout<<pattern.snp1<<":pattern.snp1 "<<pattern.snp2<<":pattern.snp2"<<endl;
    create_contingency_table_pattern2(pattern.snp1,pattern.snp2,adr_contingence2,genos,phenos_m);
    score=g_test_2SNP(contingence2);
  }
  else{
    contingence3SNP contingence3;
    contingence3SNP* adr_contingence3 = &contingence3;
    create_contingency_table_pattern3(pattern.snp1,pattern.snp2,pattern.snp3,adr_contingence3,genos,phenos_m);
    score=g_test_3SNP(contingence3);
  }
  return(score);
}*/
/*
void update(patternscore s_opt, vector<patternscore>* adr_elite_sols){
  double min_score=99999999;
  int min_elite=0;
  for(int i=0;i<(*adr_elite_sols).size();i++){
    if ((*adr_elite_sols)[i].pattern1==s_opt.pattern1 && (*adr_elite_sols)[i].pattern2==s_opt.pattern2){
      return;
    }
    if ((*adr_elite_sols)[i].pval<min_score){
      min_score=(*adr_elite_sols)[i].pval;
      min_elite=i;
    }
  }
  (*adr_elite_sols)[min_elite]=s_opt;
}*/

bool compareByLength(const patternscore &a, const patternscore &b){
  return a.pval < b.pval;
}

  vector<patternscore> identify_best_solutions(vector<patternscore> pop, int k, int n){
    vector<patternscore> best_solutions;
    std::sort(pop.begin(), pop.end(), compareByLength);
    int i=0;
    while(i<k){
      /*if (best_solutions.pval==-1)
        {continue;}
      else
        {*/best_solutions.push_back(pop[i]);
          i=i+1;}
    return(best_solutions);
    }




void update_population(vector<patternscore> children_parents, vector<patternscore>* adr_pop,int n){
      for (int i=0;i<children_parents.size();i++){
        for (int j=0;j<(*adr_pop).size();j++){
          if(children_parents[i].idparent==(*adr_pop)[j].idparent){
        //cout<<"snp"<<(*it).snp1<<(*it).snp2<<endl;
        //cout<<"size:"<<(*itp).snp1<<(*itp).snp2<<endl;
        //cout<<(*itp).idparent<<" <-parent "<<(*itp).idparent<<endl;
        (*adr_pop)[j]=children_parents[i];
        break;
        }
      }
    }
  }

/*
  patternscore hill_climbing_lc2(patternscore s_closest_neighbour, vector<patternscore> patternscoreList){
    vector<patternscore> s_neighbours = neighbours(s_closest_neighbour,patternscoreList);
    patternscore actual_s = s_closest_neighbour;
    for (int i=0;i<patternscoreList.size();i++){
      if (s_neighbours[i].pval<actual_s.pval){
        actual_s=s_neighbours[i];
      }
    }
    return(actual_s);
  }
*/

#include "./include/parameters_file_parsing.hpp"

#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

parameters_file_parsing::parameters_file_parsing()
{
    ifstream file("./PARAMETERS_PROJECT.txt");
    if(file)
    {
        string line;
        while (!file.eof())
        {
            getline(file, line);
            if (line.length() != 0 && line[0] != '#')
            {
                import_line(line);
            }
        }
    }
    else
    {
        std::cerr << "Error while opening PARAMETERS_PROJECT.txt !\n";
    }
}


////////////////////////////////////////////////////////////



void parameters_file_parsing::import_line(string const& line)
{
    vector<string> token = this->split(line, ' ');
    string const& key = token[0];
    string & value = token[1];

    if(key == "header")
        header = atoi(value.c_str());

    else if(key == "separator")
    {
        if(value == "\t")
            separator = '\t';
        else
            separator = value.at(0);
    }
    else if(key == "k")

        k = atof(value.c_str());

    else if(key == "relinking_local_search_id")
        relinking_local_search_id = value;

    else if(key == "n")
        n = atoi(value.c_str());

    else if(key == "memetique_local_search_id")
        memetique_local_search_id = value;

    else if(key == "n_it")
        n_it = atoi(value.c_str());

    else if(key == "n_pairs_selected_parents")
        n_pairs_selected_parents = atoi(value.c_str());

    else if(key == "prob_mutation")
        prob_mutation = atoi(value.c_str());

    else if(key == "best_k")
        best_k = atoi(value.c_str());


            else {}

        }
//=================================================
// parameters_file_parsing : split
//=================================================
vector<string> parameters_file_parsing::split(string const& s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim))
        tokens.push_back(item);
    return tokens;
}

//=================================================
// parameters_file_parsing : list_parameters
//=================================================
void parameters_file_parsing::list_parameters() const
{
    cout << "########### PARAMETERS ###########\n" << "header => " << header << endl
    << "separator => " << separator << endl
    << "k => " << k << endl
    << "relinking_local_search_id => " << relinking_local_search_id << endl
    << "n => " << n << endl
    << "memetique_local_search_id => " << memetique_local_search_id << endl
    << "n_it => " << n_it << endl
    << "n_pairs_selected_parents => " << n_pairs_selected_parents << endl
    << "prob_mutation => " << prob_mutation << endl
    << "best_k => " << best_k << endl
    << "#################################" << endl;
}
//=================================================
// parameters_file_parsing : update_subset_size_large
//=================================================
/*void parameters_file_parsing::update_subset_size_large(unsigned const& n_genos)
{
    if(aco_set_size == 0)
        aco_set_size = sqrt(n_genos);
}*/

/*
int main(){
 parameters_file_parsing a;
 a.list_parameters();
} */

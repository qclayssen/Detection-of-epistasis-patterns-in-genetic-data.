#include "Parameters_file_parsing.hpp"

#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

Parameters_file_parsing::Parameters_file_parsing()
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



void Parameters_file_parsing::import_line(string const& line)
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
        relinking_local_search_id = atof(value.c_str());

    else if(key == "n")
        n = value;

    else if(key == "memetique_local_search_id")
        n_smmb_aco_runs = atoi(value.c_str());

    else if(key == "n_it")
        aco_n_ants = atoi(value.c_str());

    else if(key == "n_pairs_selected_parents")
        n_pairs_selected_parents = atoi(value.c_str());

    else if(key == "prob_mutation")
        prob_mutation = atoi(value.c_str());

    else if(key == "best_k")
        best_k = atoi(value.c_str());



    else {}

}

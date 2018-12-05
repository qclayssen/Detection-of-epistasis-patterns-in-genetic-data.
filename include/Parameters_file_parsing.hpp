#ifndef PARAMETERS_FILE_PARSING_H
#define PARAMETERS_FILE_PARSING_H


#include <string>
#include <vector>


using namespace std;

class parameters_file_parsing
{
public:
    parameters_file_parsing();
    void import_line(string const& line);
    void list_parameters() const;

// Parameters given in the OPTIONS.txt file
// Reachable from any class that include the current header (Option_file_parsing.hpp)
    int header;
    char separator;
    int k;
    string relinking_local_search_id;
    unsigned n;
    string memetique_local_search_id;
    int n_it;
    unsigned n_pairs_selected_parents;
    unsigned prob_mutation;
    unsigned best_k;


private:
    vector<string> split(string const& s, char delim);

};

#endif // PARAMETERS_FILE_PARSING_H

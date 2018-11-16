#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>

#include "CSVParser.hpp"
#include "common.h"

using namespace std;
using namespace std::chrono;

int main()
{
    // Arguments
    string genos_file;
    cout<<"Genotype file path : "<<endl;
    cin>>genos_file;
    string phenos_file;
    cout<<"Phenotype file path : "<<endl;
    cin>> phenos_file;
    int header = 1;
    char separator = ',';

//  DATA IMPORTATION
    CSVParser<int> genos_csv(genos_file, separator, header);
    CSVParser<int> phenos_csv(phenos_file, separator, header);
    blas_matrix genos = genos_csv.data();
    blas_matrix phenos_m = phenos_csv.data();
    blas_column phenos(phenos_m, 0);

    cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;
    return 0;
}

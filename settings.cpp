#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 

#include "settings.h"

using namespace std;

void read(int& cas, double& xmin, double& xmax, double& ymin, double& ymax, double& Tf, int& Nx, int& Ny) {
    ifstream inputFile("data.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Unable to open input file." << endl;
        exit(1);
    }

    inputFile >> cas;
    inputFile >> xmin;
    inputFile >> xmax;
    inputFile >> ymin;
    inputFile >> ymax;
    inputFile >> Tf;
    inputFile >> Nx;
    inputFile >> Ny;
 

    inputFile.close();
}




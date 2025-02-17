#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 

#include "settings.h"

using namespace std;

<<<<<<< HEAD
void read(int& cas, double& xmin, double& xmax, double& ymin, double& ymax, double& Tf, int& Nx, int& Ny, double& alpha_robin, double& beta_robin) {
=======
void read(int& cas, double& xmin, double& xmax, double& ymin, double& ymax, double& Tf, int& Nx, int& Ny, double& alpha_robin, double& beta_robin, int& nr) {
>>>>>>> c73c8bc (Ajout du dossier Calcul_Parallele_final)
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
    inputFile >> alpha_robin;
    inputFile >> beta_robin;
<<<<<<< HEAD
=======
    inputFile >> nr;
>>>>>>> c73c8bc (Ajout du dossier Calcul_Parallele_final)
 

    inputFile.close();
}




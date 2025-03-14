#include <iostream>
#include <fstream>
#include <string>
#include "settings.h"



int main() {

    int Nx ,kmax;
    int Ny ;
    double dx;
    double dy;

    double dt;

    double x, y, xmin =0., xmax=1., CFL=0.9, ymin=0., ymax=1.;
    int cas, nr;
    double t, Tf=1. , alpha_robin , beta_robin;

    int tag;

    read( cas,  xmin,  xmax,  ymin,  ymax,  Tf,  Nx,  Ny,  alpha_robin,  beta_robin,  nr) ;

    dx = (xmax - xmin) / (Nx+1);
    dy = (ymax - ymin) / (Ny+1);


    dt = 1e-2;
    kmax=(int(Tf/dt));


    int max_i ;
    std::cout << "Entrez le nombre de proc utilisé : ";
    std::cin >> max_i;


    for (int k =0 ;k<kmax ; ++k){

    

    // Nom du fichier de sortie
    std::string outputFileName =  "Sol." + std::to_string(k+1) + ".dat";

    // Ouvre le fichier de sortie en mode d'écriture
    std::ofstream outputFile(outputFileName, std::ios::out | std::ios::binary);

    if (!outputFile.is_open()) {
        std::cerr << "Impossible d'ouvrir le fichier de sortie." << std::endl;
        return 1;
    }


    for (int i = 0; i <= max_i-1; ++i) {
        // Nom du fichier courant
        std::string fileName = "Sol." + std::to_string(cas) + "." + std::to_string(k+1) + "." + std::to_string(i) +".txt";

        // Ouvre le fichier courant en mode binaire
        std::ifstream inputFile(fileName, std::ios::binary);

        if (!inputFile.is_open()) {
            // std::cerr << "Impossible d'ouvrir le fichier " << fileName << std::endl;
            break;  
        }

        // Lit le contenu du fichier courant et écrit dans le fichier de sortie
        outputFile << inputFile.rdbuf();


        inputFile.close();
    }

    outputFile.close();

    }

    std::cout << "Concaténation des fichiers terminée." << std::endl;

    return 0;
}

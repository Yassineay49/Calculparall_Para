#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <fstream>
#include <mpi.h>

#include "Matrix.h"
#include "charge.h"
#include "functions.h"
#include "settings.h"
#include "solver.h"

using namespace std;


int main(int argc, char** argv){

    int Nx;
    int Ny;
    double dx, dy, Tf, xmin, ymin, dt;
    double t=0., xmax , ymax;
    int cas;
    double alpha_robin, beta_robin, C ;

    read(cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny,alpha_robin, beta_robin);

    std::vector<double> u0_total((Nx-1)*(Ny-1), 10.0); // Initialisation de u 
    std::vector<double> Uf((Nx-1)*(Ny-1), 0.0);
    std::vector<double> B((Nx-1)*(Ny-1), 0.0);
    std::vector<double> V((Nx-1)*(Ny-1), 0.0);


    dx = (xmax-xmin)/Nx ;
    dy = (ymax-ymin)/Ny ;
    C = ( 2 * beta_robin * dy) / alpha_robin ;

    printf("%d\n",cas); 


    // Initialisation MPI
    MPI_Init(&argc, &argv);
    int Me, Np;
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);

    // Calcul des indices de début et de fin pour chaque processus MPI
    int iBeg, iEnd;
    int val = (Nx + 1) * (Ny + 1);
    charge_a(Me, val, Np, iBeg, iEnd);


    MPI_Status status;
    int size_loc = iEnd - iBeg + 1;   



    int j=0 , Nmax=40, pas, lig;
    double eProd_scailon = 1e-4 ;
    double x,y;
    dt = 1e-5 ; 
    Tf = 200 * dt ;
    t=0. ;

    std::vector<double> vec1((Nx-1)*(size_loc-1), 0.0);
    std::vector<double> vec2((Nx-1)*(size_loc-1), 0.0);
    std::vector<double> vec3((Nx-1)*(size_loc-1), 0.0);


    while (t < Tf) {


        std::string nomFichier = "sol." + std::to_string(j) + ".dat";
        std::ofstream fichier(nomFichier);

        B = Source_term(Nx, Ny, dx, dy,  xmin,  ymin, xmax,  ymax,  t,  dt, cas, C, vec1, vec2, vec3, Me, Np);
        x= xmin +dx ;
        y=ymin +dy;
        pas = 0 ; 

        for (int i = 0; i < (Nx - 1) * (Ny - 1); ++i) {
           
            V[i] = u0_total[i] + dt * f(x,y,t,cas,xmax,ymax) - B[i];
            pas += 1;
            if (pas == Nx - 1) {
                x = xmin ;
                y = y + dy;
                pas = 0;
            }
            
            x = x + dx;
        
        }
        Uf = BiCGstab( Nx,  Ny,  dx,  dy,  xmin,  xmax,  dt,  V,  eProd_scailon,  Nmax, C, Me, Np);

        j = j + 1;

        x= xmin +dx ;
        y=ymin +dy;
        pas = 0 ; 

        for (int i = 0; i < (Nx - 1) * (Ny - 1); ++i) {
            
            fichier << x << " " << y << " " << u0_total[i] << endl;
            // fichier << x << " " << y << " " << sol_exacte(x,y,t,cas) << endl;
            lig += 1;
            pas += 1;
            if (pas == Nx - 1) {
                x = xmin ;
                y = y + dy;
                pas = 0;
            }
            
            x = x + dx;
        }
  
        t = t + dt;
        
        for (int i = 0; i < (Nx - 1) * (Ny - 1); ++i) {
            u0_total[i] = Uf[i] ;
        }
    }
   
    return 0; 
}
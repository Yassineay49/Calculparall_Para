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

    read(cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny);

    std::vector<double> u0_total((Nx-1)*(Ny-1), 10.0); // Initialisation de u 
    std::vector<double> Uf((Nx-1)*(Ny-1), 0.0);
    std::vector<double> B((Nx-1)*(Ny-1), 0.0);
    std::vector<double> V((Nx-1)*(Ny-1), 0.0);



    read(cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny);

    dx = (xmax-xmin)/Nx ;
    dy = (ymax-ymin)/Ny ;

    printf("%d\n",cas);    

    // Uf = produitmatvect( Nx,  Ny,  dx,  dy,  xmin,  ymin,  dt,u0_total );
    // B = Source_term(Nx, Ny, dx, dy,  xmin,  ymin, xmax,  ymax,  t,  cas);

    int s=0;
    for (int j = 0; j < Ny - 1; ++j) {
        for (int i = 0; i < Nx - 1; ++i) {
            
            // printf("%f\n",Uf[s]);
            // printf("%f\n",B[s]) ;

            s+=1;
        }
    }
    


    int j=0 , Nmax=40, pas, lig;
    double eProd_scailon = 1e-4 ;
    double x,y;
    dt = 1e-5 ; 
    Tf = 200 * dt ;
    t=0. ;

    while (t < Tf) {


        std::string nomFichier = "sol." + std::to_string(j) + ".dat";
        std::ofstream fichier(nomFichier);

        B = Source_term(Nx, Ny, dx, dy,  xmin,  ymin, xmax,  ymax,  t,  dt, cas);
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
        Uf = BiCGstab( Nx,  Ny,  dx,  dy,  xmin,  xmax,  dt,  V,  eProd_scailon,  Nmax );

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
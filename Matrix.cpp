#include <iostream>
#include <vector>
#include <algorithm> 
#include<fstream>
#include<mpi.h>

#include "Matrix.h"

#include "functions.h"

#include "settings.h"


using namespace std;


std::vector<double> produitmatvect(int& Nx, int& Ny, double& dx, double& dy, double& xmin, double& ymin, double& dt,std::vector<double>& u, double& C, int& Me , int& Np ) {


        std::vector<double> produit((Nx-1) * (Ny-1));

        int num_elements = (Nx-1) * (Ny-1);



        if (Me > 0 and Me < Np - 1 ){
 

            for (int i = 0; i < num_elements; i++) {

                    double alpha = 1 + (2*dt)/(dx*dx) + (2*dt/dy*dy) ;
                    double beta = -dt/(2*dx*dx) ;
                    double gamma = -dt/(2*dy*dy) ; 

                    if ( i < Nx-1)
                    {
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma) *u[i] + beta*u[i+1] + 2*gamma * u[i+Nx-1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1] + 2*gamma * u[i+Nx-1] ;

                        }
                        else{
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1] + beta*u[i+1] + 2*gamma * u[i+Nx-1];

                        }
                    }

                    if (i >= num_elements - Nx + 1)
                    {
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i+1] + 2*gamma* u[i-Nx +1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1]  + 2*gamma* u[i-Nx +1] ;

                        }
                        else{
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1] + beta*u[i+1]  + 2*gamma* u[i-Nx +1] ;

                        }
                    }


                    if (i >= Nx-1 and i < num_elements - Nx + 1 ){
                        
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha)*u[i] + beta*u[i+1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha )*u[i] + beta*u[i-1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1];

                        }
                        else{
                            produit[i] = (alpha)*u[i] + beta*u[i-1] + beta*u[i+1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1];

                        }
                    }
                    
                    
                
            }

        return produit ;

    }


    if (Me == 0 ){
 

            for (int i = 0; i < num_elements; i++) {

                    double alpha = 1 + (2*dt)/(dx*dx) + (2*dt/dy*dy) ;
                    double beta = -dt/(2*dx*dx) ;
                    double gamma = -dt/(2*dy*dy) ; 

                    if ( i < Nx-1)
                    {
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha ) *u[i] + beta*u[i+1] + gamma * u[i+Nx-1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha )*u[i] + beta*u[i-1] + gamma * u[i+Nx-1] ;

                        }
                        else{
                            produit[i] = (alpha )*u[i] + beta*u[i-1] + beta*u[i+1] + gamma * u[i+Nx-1];

                        }
                    }

                    if (i >= num_elements - Nx + 1)
                    {
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i+1] + 2*gamma* u[i-Nx +1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1]  + 2*gamma* u[i-Nx +1] ;

                        }
                        else{
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1] + beta*u[i+1]  + 2*gamma* u[i-Nx +1] ;

                        }
                    }


                    if (i >= Nx-1 and i < num_elements - Nx + 1 ){
                        
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha)*u[i] + beta*u[i+1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha )*u[i] + beta*u[i-1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1];

                        }
                        else{
                            produit[i] = (alpha)*u[i] + beta*u[i-1] + beta*u[i+1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1];

                        }
                    }
                    
                    
                
            }

        return produit ;

    }


    if (Me > 0 and Me < Np - 1 ){
 

            for (int i = 0; i < num_elements; i++) {

                    double alpha = 1 + (2*dt)/(dx*dx) + (2*dt/dy*dy) ;
                    double beta = -dt/(2*dx*dx) ;
                    double gamma = -dt/(2*dy*dy) ; 

                    if ( i < Nx-1)
                    {
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma) *u[i] + beta*u[i+1] + 2*gamma * u[i+Nx-1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1] + 2*gamma * u[i+Nx-1] ;

                        }
                        else{
                            produit[i] = (alpha - C * gamma)*u[i] + beta*u[i-1] + beta*u[i+1] + 2*gamma * u[i+Nx-1];

                        }
                    }

                    if (i >= num_elements - Nx + 1)
                    {
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha )*u[i] + beta*u[i+1] + gamma* u[i-Nx +1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha )*u[i] + beta*u[i-1]  + gamma* u[i-Nx +1] ;

                        }
                        else{
                            produit[i] = (alpha )*u[i] + beta*u[i-1] + beta*u[i+1]  + gamma* u[i-Nx +1] ;

                        }
                    }


                    if (i >= Nx-1 and i < num_elements - Nx + 1 ){
                        
                        if ((i) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha)*u[i] + beta*u[i+1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1] ;
                            
                        }
                        else if ((i+1) % (Nx-1) == 0) {
                            
                            produit[i] = (alpha )*u[i] + beta*u[i-1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1];

                        }
                        else{
                            produit[i] = (alpha)*u[i] + beta*u[i-1] + beta*u[i+1] + gamma * u[i+Nx-1] + gamma* u[i-Nx +1];

                        }
                    }
                    
                    
                
            }

        return produit ;

    }
}


std::vector<double> Source_term(int& Nx, int& Ny, double& dx, double& dy, double& xmin, double& ymin,double& xmax, double& ymax, double& t, double& dt,int& cas, double& C, std::vector<double>& vec1, std::vector<double>& vec2, std::vector<double>& vec3, int& Me, int& Np) {

    std::vector<double> B((Nx-1) * (Ny-1),0.0);
    int pas = 0;
    double y,x ;

    double beta = -dt/(2*dx*dx) ;
    double gamma = -dt/(2*dy*dy) ;


    for (int j = 0; j < Ny - 1; ++j) {
        y = ymin + (j+1) * dy;
        for (int i = 0; i < Nx - 1; ++i) {
            x = xmin + (i+1) * dx;
            if (i==0){
                B[pas] += beta * g(0,y,t,cas) ;
                
            }
            if (j==0){
                B[pas] +=  gamma * (vec1[i]  - vec3[i] + (C*vec2[i])) ;
            }
            if (i==Nx-2){
                B[pas] += beta * g(xmax,y,t,cas) ;
            }
            if (j==Ny-2){
                B[pas] += gamma * h(x,ymax,t,cas) ; 
                
            }
            pas = pas + 1;
        }
    }

    return B;
}

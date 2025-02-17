#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <fstream>
#include <mpi.h>

// Inclusions des headers depuis le dossier src/
#include "src/Matrix.h"
#include "src/charge.h"
#include "src/functions.h"
#include "src/settings.h"
#include "src/solver.h"

using namespace std;



// Programme principal pour résoudre l’équation de conduction instationnaire par la méthode de Schwarz additive
int main(int argc,char* argv[])
 {
    clock_t debut, fin ;
    int Nx,Ny,cas,kmax,Me,ibeg,iend,nproc;
    long clk_tck = CLOCKS_PER_SEC ; 
    debut=clock() ; 
    double Temps_exec ; 


    // Fichier pour enregistrer les erreurs
    std::ofstream errorFile("results/error_log.txt", std::ios::app); 


    // Lecture des paramètres globaux depuis un fichier de configuration
    int nr;
    double Lx,Ly, xmin , ymin , tf, alpha_robin , beta_robin ;
    read(cas, xmin, Lx, ymin, Ly, tf, Nx, Ny,alpha_robin, beta_robin , nr);

    // Définition des grilles en x et y
    std::vector<double> x(Nx + 2);
    std::vector<double> y(Ny + 2);

    double epsilon=0.00001 , D=1;
    double dt=1e-2 , dx,dy,tn, erreur_schwartz = 20., swhartz=1 ;
    int Nf , tag = 5;

    kmax=1000;
    dx=Lx/(Nx+1);
    dy=Ly/(Ny+1);

    double alpha=1+dt*D*((2/(pow(dx,2))+2/pow(dy,2)));
    double beta=(-1/pow(dx,2))*D*dt;
    double gama=(-1/pow(dy,2))*D*dt;
    
    int nb_shwartz=50;

    

    for(int i=0;i<Nx+2;i++)
    {
        x[i]=i*dx;
    }
    
    for(int i=0;i<Ny+2;i++)
    {
        y[i]=i*dy;
    }    
    
    // Initialisation de MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&Me);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    
    charge_a(Me,Ny,nproc,&ibeg,&iend);
    
    Nf=taille(Me,ibeg,iend,nr,nproc);
    
    // Initialisation des vecteurs pour les calculs
    std::vector<double> U0(Nx * Nf, 10.0) , b(Nx * Nf) , Uf(Nx * Nf);  
    std::vector<double> stencil1, stencil2; 

    std::vector<double> check(Nx, 0.0); 
    int nb=0 , ni = 0 , Np,Np_recouv;

    // Gestion des tailles des régions de recouvrement
    if(Me==0)
    {Np=iend-nr;}
    else
    { Np=iend-ibeg; }


    if(Me==0)
    {    Np_recouv=iend-ibeg+nr;}
    else{
        Np_recouv=iend-ibeg+2*nr;}

    // Allocation des stencils en fonction des conditions de Robin
    if (alpha_robin == 0.) { 
        stencil1.resize(Nx, 0.0); 
        stencil2.resize(Nx, 0.0);
    } else {
        stencil1.resize(3 * Nx, 0.0); 
        stencil2.resize(3*Nx, 0.0);
    }    
    
   
    tn = dt ;
    while(tn<tf)
    {

        if ( nproc == 1 ){
            // Résolution séquentielle
            b=somme_vecteur(U0,SOURCESEQ(D,dy,dx,dt,tn,Lx,Ly,Nx,Ny,x,y,U0,cas,ibeg,iend,nproc,Me),Nx*Ny);
            Uf=BICGstabSEQ(alpha,beta,gama, b , Nx , Nf ,Me,nproc,dt,D,dy,kmax,epsilon);
            U0=Uf;
        }
        else {
            while(nb <  nb_shwartz and swhartz > 1e-2)
            {

                if ( alpha_robin == 0.){
                    // Résolution parallèle avec les Conditions de raccord de Dirichlet
                    b=somme_vecteur(U0,SOURCE0(D,dy,dx,dt,tn,Lx,Ly,Nx,Ny,x,y,U0,cas,ibeg,iend,nproc,Me,nr,stencil1,stencil2),Nx*Nf);
                    Uf=BICGstab1(alpha,beta,gama, b , Nx , Nf , alpha_robin,beta_robin,Me,nproc,dt,D,dy,kmax,epsilon);
                    U0=Uf;
                
                    if(Me>0)
                    {
                        MPI_Send(&U0[(2*nr)*Nx],Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD);
                        MPI_Send(&U0[(2*nr-1)*Nx],Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD);
                        
                        MPI_Recv(&stencil1[0],Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&check[0],Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        erreur_schwartz=erreur(U0,check,Nx,0);
                
                    }
                    if(Me<nproc-1)
                    {
                        MPI_Send(&U0[Np*Nx],Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD);
                        MPI_Send(&U0[(Np+1)*Nx],Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD);
                
                        MPI_Recv(&stencil2[0],Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&check[0],Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        erreur_schwartz=erreur(U0,check,Nx,Np_recouv);
                
                    }
                }

                else {
                    // Résolution parallèle avec les Conditions de raccord de Dirichlet
                    b=somme_vecteur(U0,SOURCE1(D,dy,dx,dt,tn,Lx,Ly,Nx,Ny,x,y,U0,cas,ibeg,iend,nproc,Me,nr,alpha_robin,beta_robin,stencil1,stencil2),Nx*Nf);
                    Uf=BICGstab1(alpha,beta,gama, b , Nx , Nf , alpha_robin,beta_robin,Me,nproc,dt,D,dy,kmax,epsilon);
                    U0=Uf;

                    if(Me>0)
                    {
                        MPI_Send(&U0[(2*nr)*Nx],3*Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD);
                        MPI_Send(&U0[(2*nr-1)*Nx],Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD);

                        MPI_Recv(&stencil1[0],3*Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&check[0],Nx, MPI_DOUBLE,Me-1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        erreur_schwartz=erreur(U0,check,Nx,0);
                
                    }
                    if(Me<nproc-1)
                    {
                        MPI_Send(&U0[Np*Nx ],3*Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD);
                        MPI_Send(&U0[(Np+1)*Nx],Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD);
                
                        MPI_Recv(&stencil2[0],3*Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&check[0],Nx, MPI_DOUBLE,Me+1,tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        erreur_schwartz=erreur(U0,check,Nx,Np_recouv);

                    }
                }
                // Mise à jour des erreurs Schwarz
                MPI_Allreduce(&erreur_schwartz,&swhartz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

                // Décommentez cette partie du code si vous voulez tracer l'erreur de schwartz
                // if ( Me == 0){
                //     std :: cout << erreur_schwartz << std :: endl ;
                //     errorFile << nb << " " << erreur_schwartz << std::endl;
                // }
                nb += 1 ;   
        }
        }
        tn +=dt ;  
        ni += 1 ;


        sauvegarder_resultats(cas ,ni, Me, Nx, ibeg, iend, nproc, x ,y, U0 , nr);

    }
    
    // Finalisation MPI et affichage du temps d'exécution
    MPI_Finalize();

    fin=clock() ; 
    Temps_exec = (double)(fin-debut)/(double)clk_tck ; 
    std::cout << "\n\nTemps d'exécution : " << Temps_exec << " secondes\n\n";

    return 0;
     }
        
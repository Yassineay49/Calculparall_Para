<<<<<<< HEAD
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
=======
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


int main(int argc,char* argv[])
 {
    clock_t debut, fin ;
    int Nx,Ny,cas,kmax,Me,ibeg,iend,nproc;
    long clk_tck = CLOCKS_PER_SEC ; 
    debut=clock() ; 
    double Temps_exec ; 

    std::ofstream errorFile("error_log.txt", std::ios::app); 



    int nr;
    double Lx,Ly, xmin , ymin , tf, alpha_robin , beta_robin ;
    read(cas, xmin, Lx, ymin, Ly, tf, Nx, Ny,alpha_robin, beta_robin , nr);

    
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
    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&Me);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    
    charge_a(Me,Ny,nproc,&ibeg,&iend);
    
    Nf=taille(Me,ibeg,iend,nr,nproc);
    
    
    std::vector<double> U0(Nx * Nf, 10.0) , b(Nx * Nf) , Uf(Nx * Nf);  
    std::vector<double> stencil1, stencil2; 

    std::vector<double> check(Nx, 0.0); 
    int nb=0 , ni = 0 , Np,Np_recouv;


    if(Me==0)
    {Np=iend-nr;}
    else
    { Np=iend-ibeg; }


    if(Me==0)
    {    Np_recouv=iend-ibeg+nr;}
    else{
        Np_recouv=iend-ibeg+2*nr;}


    if (alpha_robin == 0.) { 
        stencil1.resize(Nx, 0.0); 
        stencil2.resize(Nx, 0.0);
    } else {
        stencil1.resize(3 * Nx, 0.0); 
        stencil2.resize(3*Nx, 0.0);
    }    
    
   
    tn = dt ;
    while(tn<2*dt)
    {

        if ( nproc == 1 ){
            b=somme_vecteur(U0,SOURCESEQ(D,dy,dx,dt,tn,Lx,Ly,Nx,Ny,x,y,U0,cas,ibeg,iend,nproc,Me),Nx*Ny);
            Uf=BICGstabSEQ(alpha,beta,gama, b , Nx , Nf ,Me,nproc,dt,D,dy,kmax,epsilon);
            U0=Uf;
        }
        else {
            while(nb <  nb_shwartz and swhartz > 1e-2)
            {

                if ( alpha_robin == 0.){
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
                MPI_Allreduce(&erreur_schwartz,&swhartz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                if ( Me == 0){
                    std :: cout << erreur_schwartz << std :: endl ;
                    errorFile << nb << " " << erreur_schwartz << std::endl;
                }
                nb += 1 ;   
        }
        }
        tn +=dt ;  
        ni += 1 ;


        sauvegarder_resultats(cas ,ni, Me, Nx, ibeg, iend, nproc, x ,y, U0 , nr);

    }
    
    MPI_Finalize();

    fin=clock() ; 
    Temps_exec = (double)(fin-debut)/(double)clk_tck ; 
    std::cout << "\n\nTemps d'exécution : " << Temps_exec << " secondes\n\n";

    return 0;
     }
        
>>>>>>> c73c8bc (Ajout du dossier Calcul_Parallele_final)

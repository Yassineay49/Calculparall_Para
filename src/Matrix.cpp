#include <iostream>
#include <vector>
#include <algorithm> 
#include<fstream>
#include<mpi.h>

#include "Matrix.h"

#include "functions.h"

#include "settings.h"


using namespace std;


// Fonction pour déterminer la taille de la portion de données pour chaque processus MPI
int size_proc(int rank,int ibeg,int iend,int nr,int  nproc)
{

    if(rank==0)
    {
        return (iend+nr-ibeg+1);
    }
    
    else if(rank==nproc-1)
    {
        return (iend+nr-ibeg+1);
    }

    else
    {
       return (iend-ibeg+2*nr+1);
    }
}



// Fonction pour effectuer le produit matrice-vecteur pour Np == 1
std::vector<double> produitmatvect(double a,double b,double c,std::vector<double>& x,int n,int m)
{
     std::vector<double> y(n * m);

     for(int i=0;i<n*m;i++)
     {
         y[i]=a*x[i];

         if(i%n==0)
         {
             y[i]=y[i]+b*x[i+1];

         }
         else if(i%n==n-1)
         {
             y[i]=y[i]+b*x[i-1];

         }
         else
         {
             y[i]=y[i]+b*(x[i+1]+x[i-1]);


         }
         if(i>=n)
         {
             y[i]=y[i]+c*x[i-n];

         }
         if(i<n*m-n)
         {
             y[i]=y[i]+c*x[i+n];
         }
     }


    return y;
}

// Fonction pour effectuer le produit matrice-vecteur pour Np > 1
std::vector<double> produitmatvect1(double a,double b,double c,std::vector<double>& x,int n,int m,double dy, double D, double dt, double alpha,double beta,int rank,int nproc )
{


    std::vector<double> y(n * m);


    // Produit matrice-vecteur pour les Conditions de raccord de Dirichlet
    if (alpha == 0){

     for(int i=0;i<n*m;i++)
     {
         y[i]=a*x[i];

         if(i%n==0)
         {
             y[i]=y[i]+b*x[i+1];

         }
         else if(i%n==n-1)
         {
             y[i]=y[i]+b*x[i-1];

         }
         else
         {
             y[i]=y[i]+b*(x[i+1]+x[i-1]);


         }
         if(i>=n)
         {
             y[i]=y[i]+c*x[i-n];

         }
         if(i<n*m-n)
         {
             y[i]=y[i]+c*x[i+n];
         }
     }


    return y;
    }

    // Produit matrice-vecteur pour les Conditions de raccord de Robin

    else {
        double coeff1=a+((2*beta*D*dt)/(alpha*dy));

        double coeff2=2*c;
        
        if(rank==0)
        {

        for(int i=0;i<n*m-n;i++)
        {
            y[i]=a*x[i];

            if(i%n==0)
            {
                y[i]=y[i]+b*x[i+1];

            }
            else if(i%n==n-1)
            {
                y[i]=y[i]+b*x[i-1];

            }
            else
            {
                y[i]=y[i]+b*(x[i+1]+x[i-1]);


            }
            if(i>=n)
            {
                y[i]=y[i]+c*x[i-n];

            }
            if(i<n*m-n)
            {
                y[i]=y[i]+c*x[i+n];
            }
        }
        y[n*m-n]=coeff2*x[n*m-2*n]+coeff1*x[n*m-n]+b*x[n*m-n+1];

        y[n*m-1]=coeff2*x[n*m-1-n]+coeff1*x[n*m-1]+b*x[n*m-2];

        for(int i=n*m-n+1;i<n*m-1;i++)
        {
            y[i]=coeff2*x[i-n]+b*x[i-1]+coeff1*x[i]+b*x[i+1];
        }

    }
    else if (rank==nproc-1)
    {
        for(int i=n;i<n*m;i++)
        {
            y[i]=a*x[i];

            if(i%n==0)
            {
                y[i]=y[i]+b*x[i+1];

            }
            else if(i%n==n-1)
            {
                y[i]=y[i]+b*x[i-1];

            }
            else
            {
                y[i]=y[i]+b*(x[i+1]+x[i-1]);


            }
            if(i>=n)
            {
                y[i]=y[i]+c*x[i-n];

            }
            if(i<n*m-n)
            {
                y[i]=y[i]+c*x[i+n];
            }
        }
        y[0]=2*c*(x[n])+(a+(2*beta*D*dt/(alpha*dy)))*x[0]+b*x[1];

        y[n-1]=b*x[n-2]+(a+(2*beta*D*dt/(alpha*dy)))*x[n-1]+2*c*x[2*n-1];

        for(int i=1;i<n-1;i++)
        {
            y[i]=b*x[i-1]+(a+(2*beta*D*dt/(alpha*dy)))*x[i]+b*x[i+1]+2*c*x[i+n];
        }

    }
    else
    {
        y[0]=2*c*(x[n])+(a+(2*beta*D*dt/(alpha*dy)))*x[0]+b*x[1];
        y[n-1]=b*x[n-2]+(a+(2*beta*D*dt/(alpha*dy)))*x[n-1]+2*c*x[2*n-1];
        for(int i=1;i<n-1;i++)
        {
            y[i]=b*x[i-1]+(a+(2*beta*D*dt/(alpha*dy)))*x[i]+b*x[i+1]+2*c*x[i+n];
        }
        for(int i=n;i<n*m-n;i++)
        {
            y[i]=a*x[i];

            if(i%n==0)
            {
                y[i]=y[i]+b*x[i+1];

            }
            else if(i%n==n-1)
            {
                y[i]=y[i]+b*x[i-1];

            }
            else
            {
                y[i]=y[i]+b*(x[i+1]+x[i-1]);


            }
            if(i>=n)
            {
                y[i]=y[i]+c*x[i-n];

            }
            if(i<n*m-n)
            {
                y[i]=y[i]+c*x[i+n];
            }
        }
        y[n*m-n]=2*c*(x[n*m-2*n])+(a+(2*beta*D*dt/(alpha*dy)))*x[n*m-n]+b*x[n*m-n+1];
        for(int i=n*m-n+1;i<n*m-1;i++)
        {
            y[i]=2*c*x[i-n]+b*x[i-1]+(a+(2*beta*D*dt/(alpha*dy)))*x[i]+b*x[i+1];
        }
        y[n*m-1]=2*c*(x[n*m-1-n])+(a+(2*beta*D*dt/(alpha*dy)))*x[n*m-1]+b*x[n*m-2];
 
    }
    return y;
    }

}




int taille(int rank,int ibeg,int iend,int nr,int nproc)
{

    if(rank==0)
    {
        return (iend+nr-ibeg+1);
    }
    
    else if(rank==nproc-1)
    {
        return (iend+nr-ibeg+1);
    }

    else
    {
       return (iend-ibeg+2*nr+1);
    }
}


// Source terme pour les Conditions de raccord de Dirichlet
std::vector<double> SOURCE0(double D,double dy,double dx,double dt,double t,double Lx,double Ly,int Nx,int Ny,std::vector<double>& x,std::vector<double>& y,std::vector<double>& u ,int cas,int ibeg,int iend,int nproc,int rank,int nr,std::vector<double>& stencil1,std::vector<double>& stencil2)
{
 
    int i,j,I;
      
        int N;
        N=taille(rank,ibeg,iend,nr,nproc);
        std::vector<double> b(Nx * N);

        if(rank==0)
        {
        for(j=1;j<=iend+nr+1;++j){
            if (j==1){
                for (i=1;i<=Nx;++i){
                    I = (j-1)*Nx + i-1;
                    b[I] = dt*(f(x[i],y[j],t,Lx,Ly,cas)  + D/(dy*dy)*g(x[i],0,t,Lx,Ly,cas));
                    if (i==1)
                        b[I] +=  dt*(D/(dx*dx)* h(0,dy,t,Lx,Ly,cas));
                    if (i==Nx)
                        b[I] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),dy,Lx,Ly,t,cas)); 
            }
            }
            if (j>=2 && j<iend+nr+1){
                for (i=1;i<=Nx;++i){
                    I = (j-1)*Nx + i-1;
                    b[I] =dt*f(i*dx,j*dy,t,Lx,Ly,cas);
                    if (i==1)
                    b[I] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                    if (i==Nx)
                    b[I] += dt*(D/(dx*dx) *h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }

        if (j==iend+nr+1){
            for (i=1;i<=Nx; ++i){
                I= (j-1)*Nx + i-1;
                b[I] = dt*(f(i*dx,j*dy,t,Lx,Ly,cas) + (D/(dy*dy) )*stencil2[i-1]);
                if (i==1)
                    b[I] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                if (i==Nx)
                    b[I] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }
    }
  

    return b;
        }
    
    else if(rank==nproc-1)
    {
       
        
        for(j=ibeg-nr+1;j<=iend+1;++j){
            if (j==ibeg-nr+1){
                for (i=1;i<=Nx;++i){
                    I = (j-ibeg+nr-1)*Nx + i-1;
                    b[I] = dt*(f(x[i],y[j],t,Lx,Ly,cas)  + D/(dy*dy)*stencil1[i-1]);
                    if (i==1)
                        b[I] +=  dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                    if (i==Nx)
                        b[I] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,Lx,Ly,t,cas)); 
            }
            }
            if (j>=ibeg-nr+2 && j<iend+1){
                for (i=1;i<=Nx;++i){
                    I = (j-ibeg+nr-1)*Nx + i-1;
                    b[I] =dt*f(i*dx,j*dy,t,Lx,Ly,cas);
                    if (i==1)
                    {
                    b[I] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                    }
                    if (i==Nx)
                    {
                    b[I] += dt*(D/(dx*dx) *h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
                    }
            }
        }
        if (j==iend+1){
            for (i=1;i<=Nx; ++i){
                I= (j-ibeg+nr-1)*Nx + i-1;
                b[I] = dt*(f(i*dx,j*dy,t,Lx,Ly,cas) + D/(dy*dy) * g(i*dx,(Ny+1)*dy,t,Lx,Ly,cas));

                if (i==1)
                    b[I] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                if (i==Nx)
                    b[I] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }
    }

    return b;

    }
    else 
    {

        for(j=ibeg-nr+1;j<=iend+nr+1;++j){
            if (j==ibeg-nr+1){
                for (i=1;i<=Nx;++i){
                    I = (j-ibeg+nr-1)*Nx + i-1;
                    b[I] = dt*(f(x[i],y[j],t,Lx,Ly,cas)  + (D/(dy*dy))*stencil1[i-1]);
                    if (i==1)
                        b[I] +=  dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                    if (i==Nx)
                        b[I] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,Lx,Ly,t,cas)); 
            }
            }
            if (j>=ibeg-nr+2 && j<iend+nr+1){
                for (i=1;i<=Nx;++i){
                    I = (j-ibeg+nr-1)*Nx + i-1;
                    b[I] =dt*f(i*dx,j*dy,t,Lx,Ly,cas);
                    if (i==1)
                    b[I] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                    if (i==Nx)
                    b[I] += dt*(D/(dx*dx) *h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }

        if (j==iend+nr+1){
            for (i=1;i<=Nx; ++i){
                I= (j-ibeg+nr-1)*Nx + i-1;
                b[I] = dt*(f(i*dx,j*dy,t,Lx,Ly,cas) + (D/(dy*dy))*stencil2[i-1]) ;

                if (i==1)
                    b[I] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));
                if (i==Nx)
                    b[I] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }
    }
     
    return b;

    }

}

// Source terme pour les Conditions de raccord de Robin
std::vector<double> SOURCE1(double D,double dy,double dx,double dt,double t,double Lx,
double Ly,int Nx,int Ny,std::vector<double>& x,std::vector<double>& y,std::vector<double>& u ,int cas,int iBeg,int iEnd,int Np,int rank,
int r,double alpha,double beta,std::vector<double>& stencil1, std::vector<double> stencil2)
{
    int N = size_proc(rank,iBeg,iEnd,r,Np);

    std::vector<double> b(Nx * N);        
    int i,j,k;

    if(rank==0)
    {
    
        for(j=1;j<=iEnd+r+1;++j){

            if (j==1)
            {

                for (i=1;i<=Nx;++i)
                {

                    k = (j-1)*Nx + i-1;

                    b[k] = dt*(f(x[i],y[j],t,Lx,Ly,cas)  + D/(dy*dy)*g(x[i],0,t,Lx,Ly,cas));

                    if (i==1)

                        b[k] +=  dt*(D/(dx*dx)* h(0,dy,t,Lx,Ly,cas));

                    if (i==Nx)

                        b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),dy,Lx,Ly,t,cas)); 
                }
            }

            if (j>=2 && j<iEnd+r+1)
            {
                for (i=1;i<=Nx;++i)
                {
                    k = (j-1)*Nx + i-1;

                    b[k] =dt*f(i*dx,j*dy,t,Lx,Ly,cas);

                    if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                    if (i==Nx)

                    b[k] += dt*(D/(dx*dx) *h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
                }
        }

        if (j==iEnd+r+1)
        {
            for (i=1;i<=Nx; ++i)
            {
                k= (j-1)*Nx + i-1;

                b[k] = dt*(f(i*dx,j*dy,t,Lx,Ly,cas))+((D*dt)/(dy*dy))*(stencil2[i-1+2*Nx]-stencil2[i-1])+((2*dt*D*beta)/(alpha*dy))*stencil2[i-1+Nx];

                if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                if (i==Nx)

                    b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }
    }
    
    return b;
    }
    
    else if(rank==Np-1)
    {

        
        for(j=iBeg-r+1;j<=iEnd+1;++j)
        {
            if (j==iBeg-r+1)
            {
                for (i=1;i<=Nx;++i)
                {
                    k = (j-iBeg+r-1)*Nx + i-1;

                    b[k] = dt*(f(x[i],y[j],t,Lx,Ly,cas))+((D*dt)/(dy*dy))*(stencil1[i-1]-stencil1[i-1+2*Nx])+((2*dt*D*beta)/(alpha*dy))*stencil1[i-1+Nx];

                    if (i==1)

                        b[k] +=  dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                    if (i==Nx)

                        b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,Lx,Ly,t,cas)); 
            }
            }

            if (j>=iBeg-r+2 && j<iEnd+1)
            {
                for (i=1;i<=Nx;++i)
                {
                    k = (j-iBeg+r-1)*Nx + i-1;

                    b[k] =dt*f(i*dx,j*dy,t,Lx,Ly,cas);

                    if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                    if (i==Nx)

                    b[k] += dt*(D/(dx*dx) *h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
                }
            }

        if (j==iEnd+1)
        {
            for (i=1;i<=Nx; ++i)
            {
                k= (j-iBeg+r-1)*Nx + i-1;

                b[k] = dt*(f(i*dx,j*dy,t,Lx,Ly,cas) + D/(dy*dy) * g(i*dx,(Ny+1)*dy,t,Lx,Ly,cas));

                if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                if (i==Nx)

                    b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }
    }
    return b;

    }
        else 
    {
      
           
        for(j=iBeg-r+1;j<=iEnd+r+1;++j)
        {
            if (j==iBeg-r+1)
            {
                for (i=1;i<=Nx;++i)
                {
                    k = (j-iBeg+r-1)*Nx + i-1;

                    b[k] = dt*(f(x[i],y[j],t,Lx,Ly,cas))+((D*dt)/(dy*dy))*(stencil1[i-1]-stencil1[i-1+2*Nx])+((2*dt*D*beta)/(alpha*dy))*stencil1[i-1+Nx];
                    

                    if (i==1)

                        b[k] +=  dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                    if (i==Nx)

                        b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,Lx,Ly,t,cas)); 
            }
            }
            if (j>=iBeg-r+2 && j<iEnd+r+1)
            {
                for (i=1;i<=Nx;++i)
                {
                    k = (j-iBeg+r-1)*Nx + i-1;

                    b[k] =dt*f(i*dx,j*dy,t,Lx,Ly,cas);

                    if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                    if (i==Nx)

                    b[k] += dt*(D/(dx*dx) *h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }

        
        if (j==iEnd+r+1)
        {
            for (i=1;i<=Nx; ++i)
            {
                k= (j-iBeg+r-1)*Nx + i-1;

                b[k] = dt*(f(i*dx,j*dy,t,Lx,Ly,cas))+((D*dt)/(dy*dy))*(stencil2[i-1+2*Nx]-stencil2[i-1])+((2*dt*D*beta)/(alpha*dy))*stencil2[i-1+Nx] ;

                if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                if (i==Nx)

                    b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }
    }

    return b;
    }
    }


// Source terme pour Np ==1
std::vector<double> SOURCESEQ(double D,double dy,double dx,double dt,double t,double Lx,
double Ly,int Nx,int Ny,std::vector<double>& x,std::vector<double>& y,std::vector<double>& u ,int cas,int iBeg,int iEnd,int Np,int rank)
{

    std::vector<double> b(Nx * Ny);        
    int i,j,k;


    
    for(j=1;j<=Ny; j++){

            if (j==1)
            {

                for (i=1;i<=Nx;++i)
                {

                    k = (j-1)*Nx + i-1;

                    b[k] = dt*(f(x[i],y[j],t,Lx,Ly,cas)  + D/(dy*dy)*g(x[i],0,t,Lx,Ly,cas));

                    if (i==1)

                        b[k] +=  dt*(D/(dx*dx)* h(0,dy,t,Lx,Ly,cas));

                    if (i==Nx)

                        b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),dy,Lx,Ly,t,cas)); 
                }
            }

            if (j>=2 && j<Ny)
            {
                for (i=1;i<=Nx;++i)
                {
                    k = (j-1)*Nx + i-1;

                    b[k] =dt*f(i*dx,j*dy,t,Lx,Ly,cas);

                    if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                    if (i==Nx)

                    b[k] += dt*(D/(dx*dx) *h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
                }
        }

        if (j==Ny)
        {
            for (i=1;i<=Nx; ++i)
            {
                k= (j-1)*Nx + i-1;

                b[k] = dt*(f(i*dx,j*dy,t,Lx,Ly,cas))+((D*dt)/(dy*dy))*(g(i*dx,(Ny+1)*dy,t,Lx,Ly,cas));

                if (i==1)

                    b[k] += dt*(D/(dx*dx)* h(0,j*dy,t,Lx,Ly,cas));

                if (i==Nx)

                    b[k] +=  dt*(D/(dx*dx) * h(dx*(Nx+1),j*dy,t,Lx,Ly,cas));
            }
        }
    }
    
    return b;
   
    }

// Produit matrice-vecteur pour Np == 1
std::vector<double> produitmatvectSEQ(double a,double b,double c,std::vector<double>& x,int n,int m,double dy, double D, double dt )
{


    std::vector<double> y(n * m);

     for(int i=0;i<n*m;i++)
     {
         y[i]=a*x[i];

         if(i%n==0)
         {
             y[i]=y[i]+b*x[i+1];

         }
         else if(i%n==n-1)
         {
             y[i]=y[i]+b*x[i-1];

         }
         else
         {
             y[i]=y[i]+b*(x[i+1]+x[i-1]);


         }
         if(i>=n)
         {
             y[i]=y[i]+c*x[i-n];

         }
         if(i<n*m-n)
         {
             y[i]=y[i]+c*x[i+n];
         }
     }

    return y;

}
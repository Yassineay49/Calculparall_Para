#ifndef Matrix_H
#define Matrix_H


using namespace std;

std::vector<double> produitmatvect1(double a,double b,double c,std::vector<double>& x,int n,int m,double dy, double D, double dt,double alpha,double beta,int rank,int nproc) ;
std::vector<double> produitmatvectSEQ(double a,double b,double c,std::vector<double>& x,int n,int m,double dy, double D, double dt) ; 


int taille(int rank,int ibeg,int iend,int nr,int  nproc);
int size_proc(int rank,int ibeg,int iend,int nr,int  nproc);

std::vector<double> SOURCE0(double D,double dy,double dx,double dt,double t,double Lx,double Ly,int Nx,int Ny,std::vector<double>& x,std::vector<double>& y,std::vector<double>& u ,int cas,int ibeg,int iend,int nproc,int rank,int nr,std::vector<double>& stencil1,std::vector<double>& stencil2);
std::vector<double> SOURCE1(double D,double dy,double dx,double dt,double t,double Lx,double Ly,int Nx,int Ny,std::vector<double>& x,std::vector<double>& y,std::vector<double>& u ,int cas,int iBeg,int iEnd,int Np,int rank,int r,double alpha,double beta,std::vector<double>& stencil1, std::vector<double> stencil2) ; 

std::vector<double> SOURCESEQ(double D,double dy,double dx,double dt,double t,double Lx,double Ly,int Nx,int Ny,std::vector<double>& x,std::vector<double>& y,std::vector<double>& u ,int cas,int iBeg,int iEnd,int Np,int rank) ; 

#endif
#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
using namespace std;

<<<<<<< HEAD
double Prod_sca(const std::vector<double>& v1, const std::vector<double>& v2);
double calculateNorm(const std::vector<double>& vec);
std::vector<double> BiCGstab(int Nx, int Ny, double dx, double dy, double xmin, double ymin, double dt, std::vector<double>& b, double eProd_scailon, int Nmax, double C, int& Me, int& Np);
=======
std::vector<double> somme_vecteur(const std::vector<double>& x,const std::vector<double>& y,int n);
double produit_scalaire(const std::vector<double>& X, const std::vector<double>& Y,int N);
std::vector<double> BICGstab1(double a, double b, double c, vector<double>& u , int n, int m, double alpha, double beta, int rank, int Np, double dt, double D, double dy, int kmax, double epsilon) ;
std::vector<double> BICGstabSEQ(double alpha1, double beta1, double gama1, vector<double>& b, int n, int m, int Me, int Np, double dt, double D, double dy, int kmax, double epsilon) ;
double erreur(const std::vector<double>& u,const std::vector<double>& stencil,int Nx,int p);
void sauvegarder_resultats(int cas ,int ni, int Me, int Nx, int iBeg , int iEnd , int Np , std::vector<double> x, std::vector<double> y, std::vector<double> u , int r) ;

>>>>>>> c73c8bc (Ajout du dossier Calcul_Parallele_final)

#endif

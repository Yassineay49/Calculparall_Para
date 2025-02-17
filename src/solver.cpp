#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include <cmath> 

#include "solver.h"
#include "settings.h"
#include "Matrix.h"



// Fonction pour sauvegarder les résultats dans un fichier
void sauvegarder_resultats(int cas ,int ni, int Me, int Nx, int iBeg , int iEnd , int Np , std::vector<double> x, std::vector<double> y, std::vector<double> u , int r) {
    std::string file_name = "results/Sol." + std::to_string(cas) + "." + std::to_string(ni) + "." + std::to_string(Me) + ".txt";

    std::ofstream fichier(file_name);
    if (!fichier.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier : " << file_name << std::endl;
        return;
    }
    
    int ne = iBeg - r, ns = iEnd +r ;
    if (Me == 0){
        ne = iBeg  ;
    }
    if ( Me == Np - 1)  {
        ns = iEnd ; 
    }

    for (int i = 1; i <= Nx; i++) {
        for (int j = ne+r; j < ns + 1; j++) {
            int I = (j - ne - 1) * Nx + i - 1;
            fichier << x[i] << " " << y[j] << " " << u[I] << "\n";
        }
    }
}

// Fonction pour additionner deux vecteurs
std::vector<double> somme_vecteur(const std::vector<double>& x,const std::vector<double>& y,int n)
{
    std::vector<double> h(n);
    for(int i=0;i<n;i++){
         h[i]=x[i]+y[i];
     }
     return h;
}

// Fonction pour calculer le produit scalaire de deux vecteurs
double produit_scalaire(const std::vector<double>& X, const std::vector<double>& Y,int N)
{
    double s = 0;
    for(int i=0;i<N;i++)
    {
        s=s+X[i]*Y[i];
    }
    return s;
}

// Fonction pour calculer l'erreur quadratique entre u et un stencil
double erreur(const std::vector<double>& u,const std::vector<double>& stencil,int Nx,int p)
{
    double erreur=0;
    for(int i=1;i<=Nx;i++)
    {
        int I=(p)*Nx+(i-1);
        erreur=erreur+pow(u[I]-stencil[i-1],2);
    }
    return erreur;
}


// Solveur pour le code parallèle
std::vector<double> BICGstab1(double alpha1, double beta1, double gama1, vector<double>& b, int n, int m, double alpha_rob, double beta_rob, int Me, int Np, double dt, double D, double dy, int kmax, double epsilon) 
{   
    int n1 = n ;
    n = b.size();
    vector<double> x(n, 0.0); // Initialiser x avec des zéros
    vector<double> r(n), r_tilde(n), p(n), v(n), h(n), s(n), t(n);
    double rho0, alpha, omega, rho1, beta;

    // Calculer r0 = b - Ax0 avec x0 = 0
    vector<double> Ax_0 = produitmatvect1(alpha1, beta1, gama1, x ,n1 ,m, dy, D , dt , alpha_rob , beta_rob , Me , Np );
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ax_0[i];
    }

    r_tilde = r; // Choisir r_tilde0 = r0
    rho0 = produit_scalaire(r_tilde, r ,n);
    p = r;

    for (int k = 0; k < kmax; k++) {
        // Calculer v = A * p
        v = produitmatvect1(alpha1, beta1, gama1, p ,n1 ,m, dy, D , dt , alpha_rob , beta_rob , Me , Np );
        alpha = rho0 / produit_scalaire(r_tilde, v,n);

        // Calculer h = x + alpha * p
        for (int i = 0; i < n; ++i) {
            h[i] = x[i] + alpha * p[i];
        }

        // Calculer s = r - alpha * v
        for (int i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        // Vérifier si h est une solution assez précise
        if (sqrt(produit_scalaire(s, s,n)) < epsilon) {
            x = h;
            break; // Sortir de la boucle si la solution est assez précise
        }

        // Calculer t = A * s
        t = produitmatvect1(alpha1, beta1, gama1, s ,n1 ,m, dy, D , dt , alpha_rob , beta_rob , Me , Np );
        omega = produit_scalaire(t, s,n) / produit_scalaire(t, t,n);

        // Calculer x = h + omega * s
        for (int i = 0; i < n; ++i) {
            x[i] = h[i] + omega * s[i];
        }

        // Calculer r = s - omega * t
        for (int i = 0; i < n; ++i) {
            r[i] = s[i] - omega * t[i];
        }

        // Vérifier si x est une solution assez précise
        if (sqrt(produit_scalaire(r, r,n)) < epsilon) {
            break; // Sortir de la boucle si la solution est assez précise
        }

        rho1 = produit_scalaire(r_tilde, r,n);
        beta = (rho1 / rho0) * (alpha / omega);

        // Calculer p = r + beta * (p - omega * v)
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        rho0 = rho1;
    }

    return x;
}


// Solveur pour le code séquentiel
std::vector<double> BICGstabSEQ(double alpha1, double beta1, double gama1, vector<double>& b, int n, int m, int Me, int Np, double dt, double D, double dy, int kmax, double epsilon) 
{   
    int n1 = n ;
    n = b.size();
    vector<double> x(n, 0.0); // Initialiser x avec des zéros
    vector<double> r(n), r_tilde(n), p(n), v(n), h(n), s(n), t(n);
    double rho0, alpha, omega, rho1, beta;

    // Calculer r0 = b - Ax0 avec x0 = 0
    vector<double> Ax_0 = produitmatvectSEQ(alpha1, beta1, gama1, x ,n1 ,m, dy, D , dt );
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ax_0[i];
    }

    r_tilde = r; // Choisir r_tilde0 = r0
    rho0 = produit_scalaire(r_tilde, r ,n);
    p = r;

    for (int k = 0; k < kmax; k++) {
        // Calculer v = A * p
        v = produitmatvectSEQ(alpha1, beta1, gama1, p ,n1 ,m, dy, D , dt );
        alpha = rho0 / produit_scalaire(r_tilde, v,n);

        // Calculer h = x + alpha * p
        for (int i = 0; i < n; ++i) {
            h[i] = x[i] + alpha * p[i];
        }

        // Calculer s = r - alpha * v
        for (int i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        // Vérifier si h est une solution assez précise
        if (sqrt(produit_scalaire(s, s,n)) < epsilon) {
            x = h;
            break; // Sortir de la boucle si la solution est assez précise
        }

        // Calculer t = A * s
        t = produitmatvectSEQ(alpha1, beta1, gama1, s ,n1 ,m, dy, D , dt );
        omega = produit_scalaire(t, s,n) / produit_scalaire(t, t,n);

        // Calculer x = h + omega * s
        for (int i = 0; i < n; ++i) {
            x[i] = h[i] + omega * s[i];
        }

        // Calculer r = s - omega * t
        for (int i = 0; i < n; ++i) {
            r[i] = s[i] - omega * t[i];
        }

        // Vérifier si x est une solution assez précise
        if (sqrt(produit_scalaire(r, r,n)) < epsilon) {
            break; // Sortir de la boucle si la solution est assez précise
        }

        rho1 = produit_scalaire(r_tilde, r,n);
        beta = (rho1 / rho0) * (alpha / omega);

        // Calculer p = r + beta * (p - omega * v)
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        rho0 = rho1;
    }

    return x;
}
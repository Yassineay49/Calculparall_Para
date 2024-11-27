#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include <cmath> 

#include "solver.h"
#include "settings.h"
#include "Matrix.h"



double Prod_sca(const std::vector<double>& v1, const std::vector<double>& v2) {
    double result = 0.0;
    for (size_t i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }
    return result;
}


double calculateNorm(std::vector<double>& vec) {
    double squareSum = 0.0;
    for (double x : vec) {
        squareSum += x*x ;
    }

    return std::sqrt(squareSum);
}



//Voici Bicg pour le code séquentiel 
vector<double> BiCGstab(int Nx, int Ny, double dx, double dy, double xmin, double ymin, double dt, vector<double>& b, double eProd_scailon, int Nmax )
{
    int n = b.size();
    vector<double> x(n, 0.0); // Initialiser x avec des zéros
    vector<double> r(n), r_tilde(n), p(n), v(n), h(n), s(n), t(n);
    double rho0, alpha, omega, rho1, beta;

    // Calculer r0 = b - Ax0 avec x0 = 0
    vector<double> Ax_0 = produitmatvect(Nx, Ny, dx, dy, xmin, ymin, dt, x);
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ax_0[i];
    }

    r_tilde = r; // Choisir r_tilde0 = r0
    rho0 = Prod_sca(r_tilde, r);
    p = r;

    for (int k = 0; k < Nmax; k++) {
        // Calculer v = A * p
        v = produitmatvect(Nx, Ny, dx, dy, xmin, ymin, dt, p);
        alpha = rho0 / Prod_sca(r_tilde, v);

        // Calculer h = x + alpha * p
        for (int i = 0; i < n; ++i) {
            h[i] = x[i] + alpha * p[i];
        }

        // Calculer s = r - alpha * v
        for (int i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        // Vérifier si h est une solution assez précise
        if (sqrt(Prod_sca(s, s)) < eProd_scailon) {
            x = h;
            break; // Sortir de la boucle si la solution est assez précise
        }

        // Calculer t = A * s
        t = produitmatvect(Nx, Ny, dx, dy, xmin, ymin, dt, s);
        omega = Prod_sca(t, s) / Prod_sca(t, t);

        // Calculer x = h + omega * s
        for (int i = 0; i < n; ++i) {
            x[i] = h[i] + omega * s[i];
        }

        // Calculer r = s - omega * t
        for (int i = 0; i < n; ++i) {
            r[i] = s[i] - omega * t[i];
        }

        // Vérifier si x est une solution assez précise
        if (sqrt(Prod_sca(r, r)) < eProd_scailon) {
            break; // Sortir de la boucle si la solution est assez précise
        }

        rho1 = Prod_sca(r_tilde, r);
        beta = (rho1 / rho0) * (alpha / omega);

        // Calculer p = r + beta * (p - omega * v)
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        rho0 = rho1;
    }

    return x;
}
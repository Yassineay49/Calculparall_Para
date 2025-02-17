#include <iostream>
#include <fstream>
#include <string>
#include "settings.h"
#include <cmath> 



#include "functions.h"

<<<<<<< HEAD
double f(double x, double y, double t, int cas,double xmax, double ymax){
    if (cas==1){
        return 2* ( x - x*x + y - y*y);
        // return 5. ;
    }
    if (cas==2){
        return sin(x) + cos(y);
    }
    if (cas==3){
        return exp(-(x-(xmax/2.))*(x-(xmax/2.)) ) * exp(-(y-(ymax/2.))*(y-(ymax/2.)) ) * cos(M_PI*t/2) ;
    }
}

double g(double x, double y, double t, int cas){
    if (cas==1){
        return 0.0 ; 
    }
    if (cas==2){
        return sin(x) + cos(y);
    }
    if (cas==3){
        return 0.0 ; 
    }
}

double h(double x, double y, double t, int cas){
    if (cas==1){
        return 0.0 ; 
    }
    if (cas==2){
       return sin(x) + cos(y) ; 
    }
    if (cas==3){
        return 1.0;
    }
}

double sol_exacte(double x, double y, double t, int cas){
    if (cas==1){
        return (x-x*x)*(y-y*y);
    }
}
=======

double f(double x,double y,double t,double Lx,double Ly, int cas)
{
  if (cas == 1) {
    return 2*(y-pow(y,2)+x-pow(x,2));
  }

  else if (cas == 2) {
    return sin(x)+cos(y);
  }

  else if (cas == 3) {
    double pi=4*atan(1.);
    return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos((pi*t)/2);
  }
}


double g(double x,double y,double t,double Lx,double Ly, int cas)
{
  if (cas == 1) {
    return 0;
  }

  else if (cas == 2) {
    return sin(x)+cos(y);
  }

  else if (cas == 3) {
    return 0;
  }
}


double h(double x,double y,double t,double Lx,double Ly, int cas)
{
  if (cas == 1) {
    return 0;
  }

  else if (cas == 2) {
    return sin(x)+cos(y);
  }

  else if (cas == 3) {
    return 1;
  }
}

>>>>>>> c73c8bc (Ajout du dossier Calcul_Parallele_final)

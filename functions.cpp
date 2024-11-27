#include <iostream>
#include <fstream>
#include <string>
#include "settings.h"
#include <cmath> 



#include "functions.h"

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
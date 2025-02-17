#include <iostream>
#include <fstream>
#include <string>
#include "settings.h"
#include <cmath> 



#include "functions.h"


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


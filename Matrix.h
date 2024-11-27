#ifndef Matrix_H
#define Matrix_H


using namespace std;
std::vector<double> produitmatvect(int& Nx, int& Ny, double& dx, double& dy, double& xmin, double& ymin, double& dt,std::vector<double>& u ) ;
std::vector<double> Source_term(int& Nx, int& Ny, double& dx, double& dy, double& xmin, double& ymin,double& xmax, double& ymax, double& t, double& dt, int& cas) ;


#endif
#ifndef Matrix_H
#define Matrix_H


using namespace std;
std::vector<double> produitmatvect(int& Nx, int& Ny, double& dx, double& dy, double& xmin, double& ymin, double& dt,std::vector<double>& u, double& C,int& Me, int& Np ) ;
std::vector<double> Source_term(int& Nx, int& Ny, double& dx, double& dy, double& xmin, double& ymin,double& xmax, double& ymax, double& t, double& dt, int& cas, double& C, std::vector<double>& vec1, std::vector<double>& vec2, std::vector<double>& vec3, int& Me, int& Np) ;


#endif
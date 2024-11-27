#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
using namespace std;

double Prod_sca(const std::vector<double>& v1, const std::vector<double>& v2);
double calculateNorm(const std::vector<double>& vec);
std::vector<double> BiCGstab(int Nx, int Ny, double dx, double dy, double xmin, double ymin, double dt, std::vector<double>& b, double eProd_scailon, int Nmax);

#endif

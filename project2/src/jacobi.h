#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

const double pi = 4.0*atan(1.0);

double offdiag_sq(mat& A, int n);
double norm_sq(mat& A, int n);
void get_pivot(mat& A, int n, int& k, int& l);
void rotate(mat& A, mat& V, int k, int l, int n);
void print_matrix(mat& A, int n);
void print_eigenvalues(mat& A, int n, double a, double d);

#endif
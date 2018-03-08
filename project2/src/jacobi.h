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
double get_nth_eigenvalue(mat& D, int N, int n);
int get_nth_index(mat& D, int N, int n);
void write_n_eigs(mat& D, mat& U, vec& rho, int N, int n, string filename);
void write_ground_state(mat& H, mat& H0, mat& U, mat& U0, vec& rho, int N, double omega, string filename);

#endif
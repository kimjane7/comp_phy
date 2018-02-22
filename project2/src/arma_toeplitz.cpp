// Author: Jane Kim

// This program calls armadillo functions to solve for eigenvalues 
// and eigenvectors of a tridiagonal Toeplitz matrix.

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


int main(int argc, char *argv[]){

	int N;
	double time, lambda;

	if(argc < 1){
		cout << "Input number of mesh points N." << endl;
	}

	else{
		N = atoi(argv[1]);
	}

	double h = 1.0/N;
	double hh = h*h;
	double d = 2.0/hh;
	double a = -1.0/hh;

	// set up Toeplitz matrix
	mat A = zeros<mat>(N,N);
	for(int i = 0; i < N-1; i++){
		A(i,i) = 2.0;
		A(i+1,i) = -1.0;
		A(i,i+1) = -1.0;
	}
	A(N-1,N-1) = 2.0;

	// start timer
	clock_t initial, final;
	initial = clock();

	// solve for eigenvalues and eigenvectors
	vec eigval;
	mat eigvec;
	eig_sym(eigval,eigvec,A);

	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "Computation Time = " << time << " s" << endl;

	// print to check if eigenvalues are correct
	cout << "Calculated:\tExact:" << endl;
	for(int i = 0; i < N; i++){

		lambda = d+2.0*a*cos((i+1)*pi/(N+1));

		cout << eigval(i)/hh << "\t\t" << lambda  << endl;
	}

	return 0;
}
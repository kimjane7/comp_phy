// This program uses the Jacobi 

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;

const double pi = 4.0*atan(1.0);


int main(int argc, char *argv[]){

	int N;

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

	vec eigval;
	mat eigvec;
	eig_sym(eigval,eigvec,A);


	double lambda;

	for(int i = 0; i < N; i++){

		lambda = d+2.0*a*cos((i+1)*pi/(N+1));

		cout << eigval(i)/hh << "\t" << lambda  << endl;



	}


	return 0;
}
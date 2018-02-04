//  Author: Jane Kim

//	This program approximates the solution to the 1D Poisson equation with 
//  Dirichlet boundary conditions. The algorithm assumes the matrix is tridiagonal, 
//  with 2's along the main diagonal and -1's along the two off-diagonals.

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <string> 
#include <fstream>
#include <iomanip>

using namespace std;

// test function
double f(double x) { return 100.0*exp(-10.0*x); }

// closed-form solution of -u''(x)=f(x)
double solution(double x) { return 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x); }

int main(int argc, char *argv[]){

	int max = atoi(argv[1]);	 // highest power of 10
	string filename = argv[2];   // base file name

	int n;					     // dimension of matrix
	double h, hh;		    	 // spacing, spacing squared
	double err;					 // place holder for relative error
	double factor; 				 // place holder for precalculated factor

	ofstream ofile;

	// loop through powers of 10
	for(int i = 1; i <= max; i++){

		string outfile = filename + to_string(i) + ".dat";

		n = (int) pow(10.0,i);			  
		h = 1.0/n;
		hh = h*h;

		double *x = new double[n];        // function inputs
		double *v_approx = new double[n]; // Av=b
		double *v_exact = new double[n];  // closed-form solution

		double *b = new double[n];        // Av=b
		double *d = new double[n];        // diagonal

		// set up
		for(int j = 0; j < n; j++){

			x[j] = (j+0.5)*h;
			b[j] = hh*f(x[j]);
			v_exact[j] = solution(x[j]);
			d[j]=(j+2.0)/(j+1.0);
		}

		// forward substitution
		for(int j = 1; j < n; j++){

			b[j] = b[j]+b[j-1]/d[j-1];
		}

		// backward substitution
		v_approx[n-1] = b[n-1]/d[n-1];
		for(int j = n-2; j >= 0; j--){

			v_approx[j] = (b[j]+v_approx[j+1])/d[j];
		}

		// print 
		ofile.open(outfile);
		for(int j = 0; j < n; j++){
			err = fabs((v_approx[j]-v_exact[j])/v_exact[j]);
			ofile << setw(15) << setprecision(8) << x[j];
			ofile << setw(15) << setprecision(8) << v_approx[j];
			ofile << setw(15) << setprecision(8) << v_exact[j];
			ofile << setw(15) << setprecision(8) << log10(err) << endl;
		}

		ofile.close();

		delete [] v_approx; delete [] v_exact;
		delete [] x; delete [] b; delete [] d;

	}

	return 0;
}
//  Author: Jane Kim

//	This program approximates the solution to the 1D Poisson equation, -u''(x)=f(x), with 
//  Dirichlet boundary conditions, by converting it into a system of linear equations, Av=b.
//  The algorithm (LU decomposition) assumes the matrix A is dense.

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

// test function
double f(double x) { return 100.0*exp(-10.0*x); }

// closed-form solution of -u''(x)=f(x)
double solution(double x) { return 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x); }

int main(int argc, char *argv[]){

	int max_pow;
	string filename;

	if(argc < 3){
		cout << "Input max power of 10 and file name." << endl;
	}

	else{
		max_pow = atoi(argv[1]);    // highest power of 10
		filename = argv[2];         // base file name
	}

	int n;					        // dimension of matrix
	double h, hh;		    	    // spacing, spacing squared
	double time;		            // place holder for exact time elapsed

	ofstream ofile;

	// loop through powers of 10
	for(int i = 1; i <= max_pow; i++){

		string outfile = filename + to_string(i) + ".dat";

		n = (int) pow(10.0,i);		  
		h = 1.0/(n+1);
		hh = h*h;

		mat A = zeros<mat>(n,n);    // Av=b
		vec x(n);                   // function inputs
		vec b(n); 					// Av=b
		vec u(n); 					// exact solution
		vec err(n);					// relative error

		// set up
		for(int j = 0; j < n-1; j++){

			A(j,j) = 2.0;
			A(j+1,j) = -1.0;
			A(j,j+1) = -1.0;

			x(j) = (j+1.0)*h;
			b(j) = hh*f(x(j));
		}

		A(n-1,n-1) = 2.0;
		x(n-1) = n*h;
		b(n-1) = hh*f(x(n-1));

		// start timer
		clock_t initial, final;
		initial = clock();

		// solve Av=b for v
		vec v = solve(A,b);

		// end timer
		final = clock();
		time = (final-initial)/((double) CLOCKS_PER_SEC);

		// calculate relative error
		for(int j = 0; j < n; j++){
			u(j) = solution(x(j));
			err(j) = fabs((v(j)-u(j))/u(j));
		}

		// extract max & min of relative error
		double maxerr = max(err);
		double minerr = min(err);
		double avgerr = mean(err);

		// print results to file
		ofile.open(outfile);
		ofile << "# n = " << n << endl;
		ofile << "# time used for computation = " << time << " sec" << endl;
		ofile << "# max relative error = " << maxerr << endl;
		ofile << "# min relative error = " << minerr << endl;
		ofile << "# avg relative error = " << avgerr << endl;
		ofile << "# function input (x), approx solution (v), exact solution (u), relative error" << endl;
		for(int j = 0; j < n; j++){
			ofile << setw(15) << setprecision(8) << x(j);
			ofile << setw(15) << setprecision(8) << v(j);
			ofile << setw(15) << setprecision(8) << u(j);
			ofile << setw(15) << setprecision(8) << log10(err(j)) << endl;

		}

		ofile.close();
	}

	return 0;
}
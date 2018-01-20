//  Author: Jane Kim

//	This program approximates the solution to the 1D Poisson equation with 
//  Dirichlet boundary conditions. The algorithm assumes the matrix elements
//  along the main three diagonals are different. 


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <string> 
#include <fstream>
#include <iomanip>


using namespace std;

int main(int argc, char *argv[]){

	int max = atoi(argv[1]);	 // highest power of 10
	string filename = argv[2];   // base file name
	long int n;						 // dimension of matrix
	double h, hh;		    	 // spacing, spacing squared
	double f;					 // place holder for function value
	double err;					 // place holder for relative error

	ofstream ofile;

	// loop through powers of 10
	for(int i = 1; i <= max; i++){

		string outfile = filename + to_string(i) + ".dat";

		n = (int) pow(10.0,i); // is int large enough for n?
		h = 1.0/n;
		hh = h*h;

		double *x = new double[n]; // function inputs
		double *a = new double[n]; // lower diagonal
		double *b = new double[n]; // diagonal 
		double *c = new double[n]; // upper diagonal
		double *g = new double[n]; // Av=g
		double *v = new double[n]; 
		double *v_exact = new double[n];

		// set-up
		for(int j = 0; j < n; j++){

			// fill out vectors
			x[j] = (j+0.5)*h;
			f = 100.0*exp(-10.0*x[j]);
			g[j] = hh*f;
			v_exact[j] = 1.0-(1.0-exp(-10.0))*x[j]-exp(-10*x[j]);
			    
			// fill out matrix diagonals
			// read from file for general matrix?
			// WATCH OUT FOR ENDS OF OFF-DIAG
			a[j] = -1.0;
			b[j] = 2.0;
			c[j] = -1.0;

		}
		a[0]=0;
		c[n-1]=0;

		// forward substitution
		c[0] = c[0]/b[0];
		g[0] = g[0]/b[0];
		for(int j = 1; j < n; j++){
			c[j] = c[j]/(b[j]-c[j-1]*a[j]);
			g[j] = (g[j]-g[j-1]*a[j])/(b[j]-c[j-1]*a[j]);
		}

		// backward substitution
		v[n-1] = g[n-1];
		for(int j = n-2; j > 0; j--){
			v[j] = g[j]-c[j]*v[j+1];
		}

		// print 
		ofile.open(outfile);
		for(int j = 0; j < n; j++){
			err = fabs((v[j]-v_exact[j])/v_exact[j]);
			ofile << setw(15) << setprecision(8) << x[j];
			ofile << setw(15) << setprecision(8) << v[j];
			ofile << setw(15) << setprecision(8) << v_exact[j];
			ofile << setw(15) << setprecision(8) << log10(err) << endl;
		}
		ofile.close();

	}
	
	


	return 0;
}
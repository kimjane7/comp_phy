//  Author: Jane Kim

//	This program approximates the solution to the 1D Poisson equation with 
//  Dirichlet boundary conditions. This algorithm assumes the matrix is identically
//  -2 on the main diagonal and identically 1 on the off-diagonals.


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <string> 
#include <fstream>
#include <iomanip>


using namespace std;

int main(int argc, char *argv[]){

	int max = atoi(argv[1]);	 // highest power of 10
	string filename = argv[2];   // base file name
	long int n;					 // dimension of matrix
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
		double *b = new double[n]; // Av=b
		double *d = new double[n]; 
		double *v = new double[n]; 
		double *v_exact = new double[n]; 

		// set-up
		for(int j = 0; j < n; j++){

			// fill out vectors
			x[j] = (j+0.5)*h;
			f = 100.0*exp(-10.0*x[j]);
			b[j] = hh*f;
			v_exact[j] = 1.0-(1.0-exp(-10.0))*x[j]-exp(-10*x[j]);
			d[j] = (j+1.0)/j;
			cout << d[j] << endl;
		}
		d[0] = d[n-1] = 2.0;

		// forward subsitution
		for(int j = 2; j < n; j++) b[j] = b[j]+b[j-1]/d[i-1] ;

		// backward substitution
		for(int j = n-2; j > 0; j--) v[j] = (b[j]+v[j+1])/d[i] ;

		// print 
		ofile.open(outfile);
		ofile << "# x, approximate solution, exact solution, log(error)" << endl;
		for(int j = 0; j < n; j++){
			err = fabs((v[j]-v_exact[j])/v_exact[j]);
			ofile << setw(15) << setprecision(8) << d[j];
			ofile << setw(15) << setprecision(8) << v[j];
			ofile << setw(15) << setprecision(8) << v_exact[j];
			ofile << setw(15) << setprecision(8) << log10(err) << endl;
		}
		ofile.close();	
		


	}
	

	return 0;
}
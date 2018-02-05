//  Author: Jane Kim

//	This program approximates the solution to the 1D Poisson equation, -u''(x)=f(x), with 
//  Dirichlet boundary conditions, by converting it into a system of linear equations, Av=b.
//  The algorithm assumes A is a symmetric tridiagonal matrix.

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;

// test function
double f(double x) { return 100.0*exp(-10.0*x); }

// closed-form solution of -u''(x)=f(x)
double solution(double x) { return 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x); }

int main(int argc, char *argv[]){

	int max;
	string filename;

	if(argc < 3){
		cout << "Input max power of 10 and file name." << endl;
	}

	else{
		max = atoi(argv[1]);	  // highest power of 10
		filename = argv[2];       // base file name
	}

	int n;					      // dimension of matrix
	double h, hh;		    	  // spacing, spacing squared
	double time;		          // place holder for time elapsed
	double factor; 				  // place holder for precalculated factor

	ofstream ofile;

	// loop through powers of 10
	for(int i = 1; i <= max; i++){

		string outfile = filename + to_string(i) + ".dat";

		n = (int) pow(10.0,i);			  
		h = 1.0/(n+1);
		hh = h*h;

		double *x = new double[n];         // function inputs
		double *v = new double[n];         // approx solution of Av=b
		double *u = new double[n];         // closed-form solution
		double *err = new double[n];       // relative error

		double *b = new double[n];         // Av=b
		double *c = new double[n];         // lower and upper diagonals
		double *d = new double[n];         // diagonal

		// set up
		for(int j = 0; j < n; j++){

			x[j] = (j+1.0)*h;
			b[j] = hh*f(x[j]);
			u[j] = solution(x[j]);

			c[j] = -1.0;
			d[j] = 2.0;
		}

		// start timer
		clock_t initial, final;
		initial = clock();

		// forward substitution
		for(int j = 1; j < n; j++){
			factor = c[j-1]/d[j-1];
			d[j] = d[j]-factor*c[j-1];
			b[j] = b[j]-factor*b[j-1];
		}

		// backward substitution
		v[n-1] = b[n-1]/d[n-1];
		for(int j = n-2; j >= 0; j--){
			v[j] = (b[j]-c[j]*v[j+1])/d[j];
		}

		// end timer
		final = clock();
		time = (final-initial)/((double) CLOCKS_PER_SEC);

		// calculate max, min, avg relative errors
		double maxerr = -100.0;
		double minerr = 100.0;
		double avgerr = 0;
		for(int j = 0; j < n; j++){
			err[j] = fabs((v[j]-u[j])/u[j]);
			if(err[j] > maxerr) maxerr = err[j];
			if(err[j] < minerr) minerr = err[j];
			avgerr += err[j];
		}
		avgerr = avgerr/n;	

		// print results to file
		ofile.open(outfile);
		ofile << "# n = " << n << endl;
		ofile << "# time used for computation = " << time << " sec" << endl;
		ofile << "# max relative error = " << maxerr << endl;
		ofile << "# min relative error = " << minerr << endl;
		ofile << "# avg relative error = " << avgerr << endl;
		ofile << "# function input (x), approx solution (v), exact solution (u), log(relative error)" << endl;
		for(int j = 0; j < n; j++){
			ofile << setw(15) << setprecision(8) << x[j];
			ofile << setw(15) << setprecision(8) << v[j];
			ofile << setw(15) << setprecision(8) << u[j];
			ofile << setw(15) << setprecision(8) << log10(err[j]) << endl;
		}

		ofile.close();

		delete [] x; delete [] v; delete [] u; delete [] err;
		delete [] b; delete [] c; delete [] d;

	}

	return 0;
}
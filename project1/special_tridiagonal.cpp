//  Author: Jane Kim

//	This program approximates the solution to the 1D Poisson equation with 
//  Dirichlet boundary conditions. This algorithm assumes the matrix has
//  identical elements along the diagonal and identical (but different!) 
//  elements along the off-diagonals. 


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>


int main(int argc, char *argv[]){

	int max = atoi(argv[1]); // highest power of 10
	int n;					 // dimension of equation
	double h;				 // spacing
	double f;				 // place holder for function value
	double x;				 // place holder for function input


	// loop through powers of 10
	for(int i = min; i <= max; i++){

		n = (int) pow(10.0,i); // is int large enough for n?
		h = 1.0/n;
		hh = h*h;
		

		double *d = new double[n]; // diagonal of A
		double *l = new double[n]; // lower diagonal of A
		double *b = new double[n]; // Av=b
		double *v = new double[n]; 
		double *v_exact = new double[n]; 

		// set-up
		for(int j = 0; j < n; j++){

			// fill out b and v_exact vectors
			x = j*h;
			f = 100.0*exp(-10.0*x);
			b[j] = hh*f;
			v_exact[j] = 1.0-(1.0-exp(-10.0))*x-exp(-10*x);
			    
			// fill out matrix diagonals
			// read from file for general matrix?
			u[j] = -1.0;
			d[j] = 2.0;
			l[j] = -1.0;

		}

		// forward substitution
		u[0] = u[0]/d[0];
		b[0] = b[0]/d[0];
		for(int j = 1; j < n; j++){
			u[j] = u[j]/(d[j]-l[j]*u[j-1]);
			b[j] = (b[j]-l[j]*b[j-1])/(d[j]-l[j]*u[j-1]);
		}

		// backward substitution
		v[n-1] = b[n-1];
		for(int j = n-2; j > 0; j--){
			v[j] =
		}

	}
	

	return 0;
}
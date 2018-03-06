// This file contains the definitions of functions used to calculate 
// the lowest k eigenvalues of a real, symmtric, tridiagonal matrix 
// by the method of bisection.


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

// calculates q_k(x) for k = 1, ..., n-1
double q(int k, double x, vec& d, vec& a){

	double q = d(0)-x;
	for(int i = 1; i <= k; i++){
		q = d(i)-x-a(i)*a(i)/q;
	}

	return q;
}

// uses bisection method to find root of q_k(x) in [min, max]
double get_root(int k, double min, double max, vec& d, vec& a){

	double tolerance = 1.0E-5;
	double mid, qmin, qmid, qmax;

	while(fabs(max-min) > tolerance){

		mid = 0.5*(min+max);

		if(k==7) {
			cout << setw(15) << setprecision(8) << qmin;
			cout << setw(15) << setprecision(8) << qmid;
			cout << setw(15) << setprecision(8) << qmax << endl;
		}

		
	
 		qmin = q(k, min, d, a);
		qmid = q(k, mid, d, a);
		qmax = q(k, max, d, a);

		if(qmid == 0.0) return mid;
		if(qmin*qmid < 0.0) max = mid;
		if(qmid*qmax < 0.0) min = mid;
	}

	return 0.5*(min+max);
}

void get_eigenvalues(vec& x, vec& d, vec& a, double xmin, double xmax, int N){

	double x0, xi, offset = 0.2;
	x(0) = d(0);

	for(int i = 1; i < N; i++){
		cout << i << endl;
		x0 = get_root(i, xmin+offset, x(0)-offset, d, a);
		xi = get_root(i, x(i-1)+offset, xmax-offset, d, a);
		for(int j = 1; j < i; j++){
			x(j) = get_root(i, x(j-1)+offset, x(j)-offset, d, a);
		}
		x(0) = x0;
		x(i) = xi;
	}

	for(int i = 0; i < N; i++){
		cout << i << "\t" << x(i) << endl;
	}
}

int main(int argc, char *argv[]){

	int N;

	if(argc < 1){ cout << "Input number of mesh points N." << endl; }
	else{ 
		N = atoi(argv[1]);
	}	

	double rhomin = 0.0, rhomax = 10.0;
	double h = (rhomax-rhomin)/N, hh = h*h;
	vec rho(N), V(N);
	vec d(N), a(N), aa(N);
	vec eigvals = zeros<vec>(N);

	for(int i = 0; i < N; i++){
		rho(i) = rhomin + i*h;
		V(i) = rho(i)*rho(i);
		d(i) = 2.0/hh+V(i);
		a(i) = -1.0/hh;
	}
	a(0) = 0.0;

	double xmin, xmax, deltax;
	xmin = d(N-1)-fabs(a(N-1));
	xmax = d(N-1)+fabs(a(N-1));

	for(int i = N-2; i >= 0; i--){
		deltax = fabs(a(i)) + fabs(a(i+1));
		if(d(i)+deltax > xmax) xmax = d(i)+deltax;
		if(d(i)-deltax < xmin) xmin = d(i)-deltax;
	}

	//get_eigenvalues(eigvals, d, a, xmin, xmax, N);




	return 0;
}
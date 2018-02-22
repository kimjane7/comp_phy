// Author: Jane Kim

// This program solves for the eigenvalues of a real, symmetric matrix
// using the Jacobi rotation algorithm. 

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

double offdiag_sq(mat& A, int n){

	double sum = 0;

	// sum of squares of elements above the main diagonal
	for(int i = 0; i < n; i++){
		for(int j = i+1; j < n; j++){
			sum += A(i,j)*A(i,j);
		}
	}

	// A is symmetric
	return 2.0*sum;
}

double norm_sq(mat& A, int n){

	double norm_sq = 0;

	// sum of squares of main diagonal elements
	for(int i = 0; i < n; i++){
		norm_sq += A(i,i)*A(i,i);
	}

	// add squares of off-diagonal elements
	norm_sq += offdiag_sq(&A,n);

	return norm_sq;
}

// stores the index of the max off-diagonal element of each column
// only need to change p(k) and p(l) after each rotation
void get_pivots(mat &A, ivec &p, int k, int l, int n){

	double kmax, lmax, Aik, Ail;

	for(int i = 0; i < n; i++){

		Aik = fabs(A(i,k));
		Ail = fabs(A(i,l));

		if( Aik > kmax ) { kmax = Aik; }
		if( Ail > lmax ) { lmax = Ail; }
	}

	p(k) = kmax;
	p(l) = lmax;
}

void rotate(mat &A, mat &V, ivec &p, int k, int l, int n){

	// zero out A(k,l)=A(l,k)
	if( A(k,l) != 0.0 ){

		double c, s, t, tau;
		double cc, ss, cs;
		double Aik, Ail, Vik, Vil;
		double Akk = A(k,k), All = A(l,l), Akl = A(k,l);

		// calculate angle of rotation
		tau = 0.5*(All-Akk)/Akl;
		if(tau >= 0.0){ t = tau+sqrt(1.0+tau*tau); }
		else{ t = tau-sqrt(1.0+tau*tau); }

		cc = 1.0/(1.0+t*t);
		ss = 1.0-cc;
		cs = t*cc;
		c = sqrt(cc);
		s = t*c;

		// perform rotation
		A(k,l) = A(l,k) = 0;
		A(k,k) = cc*Akk+ss*All-2.0*cs*Akl;
		A(l,l) = ss*Akk+cc*All+2.0*cs*Akl;
		for(int i = 0; i < n; i++){

			if( (i!=k) && (i!=l) ){
				Aik = A(i,k);
				Ail = A(i,l);
				A(i,k) = c*Aik-s*Ail;
				A(i,l) = c*Ail+s*Aik;
				A(k,i) = Aik;
				A(l,i) = Ail;
			}

			// rotate eigenvectors 
			Vik = V(i,k);
			Vil = V(i,l);
			V(i,k) = c*Vik-s*Vil;
			V(i,l) = s*Vik+c*Vil;
		}
	}

	else( cout << "These elements are already zero!" << endl; )
}


int main(int argc, char *argv[]){

	int N;

	if(argc < 1){ cout << "Input dimension of matrix." << endl; }
	else{ N = atoi(argv[1]); }

	double h = 1.0/N, hh = h*h;
	double time, epsilon;

	// set up matrix
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

	// solve for eigenvalues
	double tau, t, s, c;



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
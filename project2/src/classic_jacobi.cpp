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
	norm_sq += offdiag_sq(A,n);

	return norm_sq;
}

void get_pivots(mat& A, ivec& p, int n){

	int jmax; 		      // placeholder for index of largest element
	double imax = 0.0;    // value of max element in row i
	double Aij;

	// p(i) = index of largest off-diag elem in row i
	// A is symmetric, so only consider elements above diagonal
	for(int i = 0; i < n-1; i++){
		for(int j = i+1; j < n; j++){
			Aij = fabs(A(i,j));
			if(Aij > imax){
				imax = Aij; 
				jmax = j;
			}
		}
		p(i) = jmax;
		imax = 0.0;
	}
}

void update_pivots(mat& A, ivec& p, int k, int l, int n){

	double kmax = 0.0, lmax = 0.0;
	double Akj, Alj;

	for(int i = 0; i < n; i++){
		if(i = k){
			for(int j = l+1; j < n; j++){
				Akj = fabs(A(k,j));
				if(Akj > kmax) kmax = j;
			}
			p(k) = kmax;
		}
		if( (i=l) && (l<n-1))
	}



	// l may equal n-1
	if(l < n-1){
		for(int j = l+1; j < n; j++){
			Alj = fabs(A(l,j));
			if(Alj > lmax) lmax = j;
		}
		p(l) = lmax;	
	}
}

void rotate(mat& A, mat& V, ivec& p, int k, int l, int n){

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
		A(k,l) = 0.0;
		A(l,k) = 0.0;
		A(k,k) = cc*Akk+ss*All-2.0*cs*Akl;
		A(l,l) = ss*Akk+cc*All+2.0*cs*Akl;
		for(int i = 0; i < n; i++){

			if( (i!=k) && (i!=l) ){
				Aik = A(i,k);
				Ail = A(i,l);
				A(i,k) = c*Aik-s*Ail;
				A(i,l) = c*Ail+s*Aik;
				A(k,i) = A(i,k);
				A(l,i) = A(i,l);
			}

			// rotate eigenvectors 
			Vik = V(i,k);
			Vil = V(i,l);
			V(i,k) = c*Vik-s*Vil;
			V(i,l) = s*Vik+c*Vil;
		}
	}

	else{ cout << "These elements are already zero!" << endl; } 
}

void print_matrix(mat& A, int n){

	for(int i = 0; i < n; i++){
		cout << "[";
		for(int j = 0; j < n; j++){
			cout << setw(10) << setprecision(3) << A(i,j);
		}
		cout << "]" << endl;
	}
	cout << "\n" << endl;
}

int main(int argc, char *argv[]){

	int n, k, l, iterations = 0;
	double max_offdiag, Aij; 

	if(argc < 1){ cout << "Input dimension of matrix." << endl; }
	else{ n = atoi(argv[1]); }

	double h = 1.0/n, hh= h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double time, lambda, epsilon = 10E-10;

	// set-up matrix to diagonalize
	mat A = zeros<mat>(n,n);
	for(int i = 0; i < n-1; i++){
		A(i,i) = d;
		A(i+1,i) = a;
		A(i,i+1) = a;
	}
	A(n-1,n-1) = d;

	print_matrix(A,n);

	// set-up matrix of eigenvectors
	mat V = eye<mat>(n,n);

	// set-up list of indices of the max off-diag elem of each column
	ivec p(n-1);
	get_pivots(A, p, n);

	// start timer
	clock_t initial, final;
	initial = clock();

	// solve for eigenvalues and eigenvectors
	while(offdiag_sq(A,n) > epsilon){

		// find indices (k,l) of max pivot
		max_offdiag = 0.0;
		for(int i = 0; i < n-1; i++){
			Aij = fabs(A(i,p(i)));
			if(Aij > max_offdiag){
				max_offdiag = Aij;
				k = i;
				l = p(i);
			}
		}

		cout << "(k,l) = " << k << "\t" << l << endl;

		// apply rotation
		rotate(A, V, p, k, l, n);

		print_matrix(A,n);

		// update list of pivots
		update_pivots(A, p, k, l, n);

		iterations += 1;

	}


	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "Diagonalized in " << iterations << " iterations" << endl;
	cout << "Computation Time = " << time << " s" << endl;

	// print to check if eigenvalues are correct
	cout << "Calculated:\tExact:" << endl;
	for(int i = 0; i < n; i++){

		lambda = d+2.0*a*cos((i+1)*pi/(n+1));

		cout << A(i,i) << "\t\t" << lambda  << endl;
	}

	return 0;
}
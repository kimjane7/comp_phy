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

void get_pivot(mat& A, int n, int& k, int& l){

	double Aij,max_offdiag = 0.0;

	for(int i = 0; i < n-1; i++){
		for(int j = i+1; j < n; j++){
			Aij = fabs(A(i,j));
			if(Aij > max_offdiag){
				max_offdiag = Aij;
				k = i;
				l = j;
			}
		}
	}
}

void rotate(mat& A, mat& V, int k, int l, int n){

	// zero out A(k,l)=A(l,k)
	if( A(k,l) != 0.0 ){

		double c, s, t, tau;
		double cc, ss, cs;
		double Aik, Ail, Vik, Vil;
		double Akk = A(k,k), All = A(l,l), Akl = A(k,l);

		// calculate angle of rotation
		tau = 0.5*(All-Akk)/Akl;
		if(tau >= 0.0){ t = -tau+sqrt(1.0+tau*tau); }
		else{ t = -tau-sqrt(1.0+tau*tau); }

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
	cout << endl;
}

void print_eigenvalues(mat& A, int n, double a, double d){

	vec eigvals(n);
	double lambda;

	for(int i = 0; i < n; i++){
		eigvals(i) = A(i,i);
	}

	sort(eigvals.begin(),eigvals.end());

	// check if eigenvalues are correct
	cout << left << "* EIGENVALUES *\n";
	cout << showpoint;
	cout << setw(15) << "Calculated:";
	cout << setw(15) << "Exact:" << endl;

	for(int i = 0; i < n; i++){

		lambda = d+2.0*a*cos((i+1)*pi/(n+1));

		cout << setprecision(10) << setw(15) << eigvals(i); 
		cout << setprecision(10) << setw(15) << lambda << endl;
	}

	cout << endl;
}

int main(int argc, char *argv[]){

	int n, k, l, iterations = 0;
	double max_offdiag, Aij; 

	if(argc < 1){ cout << "Input dimension of matrix." << endl; }
	else{ n = atoi(argv[1]); }

	double h = 1.0/n, hh= h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double time, epsilon = 1E-5;

	// set-up matrix to diagonalize
	mat A = zeros<mat>(n,n);
	for(int i = 0; i < n-1; i++){
		A(i,i) = d;
		A(i+1,i) = a;
		A(i,i+1) = a;
	}
	A(n-1,n-1) = d;

	// set-up matrix of eigenvectors
	mat V = eye<mat>(n,n);
	

	// start timer
	clock_t initial, final;
	initial = clock();

	// solve for eigenvalues and eigenvectors
	while(offdiag_sq(A,n) > epsilon){

		get_pivot(A, n, k, l);

		rotate(A, V, k, l, n);

		// print_matrix(A,n);

		iterations += 1;
	}

	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "\nComputation Time = " << time << " s\n";

	cout << "Diagonalized in " << iterations << " iterations\n" << endl;

	print_eigenvalues(A,n,a,d);

	// cout << "* EIGENVECTORS *" << endl;
	// print_matrix(V,n);


	return 0;
}
// This file contains the definitions of functions used to implement 
// the Jacobi rotation algorithm for a real, symmetric matrix.

#include "jacobi.h"

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

	else{ cout << "ERROR: These elements are already zero!" << endl; } 
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


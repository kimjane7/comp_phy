// This file contains the definitions of functions used to implement 
// the Jacobi rotation algorithm for a real, symmetric matrix.

#include "jacobi.h"

double offdiag_sq(mat& A, int N){

	double sum = 0;

	// sum of squares of elements above the main diagonal
	for(int i = 0; i < N; i++){
		for(int j = i+1; j < N; j++){
			sum += A(i,j)*A(i,j);
		}
	}

	// A is symmetric
	return 2.0*sum;
}

double norm_sq(mat& A, int N){

	double norm_sq = 0;

	// sum of squares of main diagonal elements
	for(int i = 0; i < N; i++){
		norm_sq += A(i,i)*A(i,i);
	}

	// add squares of off-diagonal elements
	norm_sq += offdiag_sq(A,N);

	return norm_sq;
}

void get_pivot(mat& A, int N, int& k, int& l){

	double Aij,max_offdiag = 0.0;

	for(int i = 0; i < N-1; i++){
		for(int j = i+1; j < N; j++){
			Aij = fabs(A(i,j));
			if(Aij > max_offdiag){
				max_offdiag = Aij;
				k = i;
				l = j;
			}
		}
	}
}

void rotate(mat& A, mat& V, int k, int l, int N){

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
		for(int i = 0; i < N; i++){

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

void jacobi(mat& A, mat& V, int N){

	double epsilon = 1E-8;
	int k, l, iterations = 0;

	while(offdiag_sq(A,N) > epsilon){

		get_pivot(A, N, k, l);

		rotate(A, V, k, l, N);

		iterations += 1;
	}

	cout << "Diagonalized in " << iterations << " iterations" << endl;

}

void print_matrix(mat& A, int N){

	for(int i = 0; i < N; i++){
		cout << "[";
		for(int j = 0; j < N; j++){
			cout << setw(10) << setprecision(3) << A(i,j);
		}
		cout << "]" << endl;
	}
	cout << endl;
}

// returns nth largest eigenvalue of diagonalized matrix
double get_nth_eigenvalue(mat& D, int N, int n){

	vec eigvals(N);

	for(int i = 0; i < N; i++){
		eigvals(i) = D(i,i);
	}
	sort(eigvals.begin(),eigvals.end());

	return eigvals(n);
}

// returns index of nth largest eigenvalue
int get_nth_index(mat& D, int N, int n){

	int index;
	vec eigvals(N);

	for(int i = 0; i < N; i++){
		eigvals(i) = D(i,i);
	}
	sort(eigvals.begin(),eigvals.end());

	for(int i = 0; i < N; i++){
		if(eigvals(n) == D(i,i)) index = i;
	}

	return index;
}

// print eigenvectors corresponding to lowest n eigenvalues
void write_n_eigs(mat& D, mat& U, vec& rho, int N, int n, string filename){


	ofstream ofile;
	double lambda, u, uu;
	vec eigvals(N);
	ivec index(n);

	// sort eigenvalues by ascending order
	for(int i = 0; i < N; i++){
		eigvals(i) = D(i,i);
	}
	sort(eigvals.begin(),eigvals.end());

	// get indices of lowest 3 eigenvalues
	for(int i = 0; i < n; i++){
		for(int j = 0; j < N; j++){
			if(eigvals(i) == D(j,j)) index(i) = j;
		}
	}

	// print u(rho) for lowest three eigenvalues
	for(int j = 0; j < n; j++){

		lambda = 3.0+j*4.0;

		string outfile = filename + to_string(j+1) + ".dat";
		cout << "Writing to '" << outfile << "'... ";

		ofile.open(outfile);
		ofile << "# calculated eigenvalue = " << eigvals(j) << endl;
		ofile << "# exact eigenvalue = " << lambda << endl;
		ofile << "# rho = r/alpha, eigenvector = u(rho)" << endl;

		for(int i = 0; i < N; i++){
			u = U(i,index(j));
			uu = u*u;
			ofile << rho(i) << "\t" << u << "\t" << uu << endl;
		}
		ofile.close();
		
		cout << "done!" << endl;
	}
}

void write_ground_state(mat& H, mat& H0, mat& U, mat& U0, vec& rho, int N, double omega, string filename){

	double u, u0;

	cout << "Writing to '" << filename << "'... ";

	ofstream ofile;
	ofile.open(filename);

	ofile << "# N = " << N << endl;
	ofile << "# omega_r = " << omega << endl;
	ofile << "# computation time = " << time << endl;
	ofile << "# first eigenvalue [Coulomb] = " << get_nth_eigenvalue(H, N, 0) << endl;
	ofile << "# first eigenvalue [no Coulomb] = " << get_nth_eigenvalue(H0, N, 0) << endl;
	ofile << "# rho, |u(rho)|^2 [Coulomb], |u0(rho)|^2 [no Coulomb]" << endl;
	
	int j = get_nth_index(H,N,0);
	int j0 = get_nth_index(H0,N,0);
	for(int i = 0; i < N; i++){
		u = U(i,j);
		u0 = U0(i,j0);
		ofile << setw(15) << setprecision(8) << rho(i);
		ofile << setw(15) << setprecision(8) << u*u;
		ofile << setw(15) << setprecision(8) << u0*u0 << endl;
	}

	ofile.close();

	cout << "done!" << endl;
}

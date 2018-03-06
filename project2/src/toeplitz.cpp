// g++ -std=c++11 toeplitz.cpp jacobi.cpp -o toeplitz -larmadillo

#include "jacobi.h"

void check_eigenvalues(mat& A, int n, double a, double d){

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

	int n;

	if(argc < 1){ cout << "Input dimension of matrix." << endl; }
	else{ n = atoi(argv[1]); }

	double h = 1.0/n, hh= h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double time, epsilon = 1E-5;
	double max_offdiag, Aij; 
	int k, l, iterations = 0;

	// set-up toeplitz matrix to diagonalize
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

		iterations += 1;
	}

	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "\nComputation Time = " << time << " s\n";

	cout << "Diagonalized in " << iterations << " iterations\n" << endl;

	check_eigenvalues(A,n,a,d);

	cout << "* EIGENVECTORS *" << endl;
	print_matrix(V,n);


	return 0;
}
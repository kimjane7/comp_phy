#include "jacobi.h"


int main(int argc, char *argv[]){

	int N;

	if(argc < 1){ cout << "Input dimension of matrix." << endl; }
	else{ N = atoi(argv[1]); }

	double h = 1.0/(N+1), hh= h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double time, epsilon = 1E-5;
	double max_offdiag, Aij; 
	int k, l, iterations = 0;

	// set-up toeplitz matrix to diagonalize
	mat A = zeros<mat>(N,N);
	for(int i = 0; i < N-1; i++){
		A(i,i) = d;
		A(i+1,i) = a;
		A(i,i+1) = a;
	}
	A(N-1,N-1) = d;

	// set-up matrix of eigenvectors
	mat V = eye<mat>(N,N);

	// start timer
	clock_t initial, final;
	initial = clock();

	// solve for eigenvalues and eigenvectors
	while(offdiag_sq(A,N) > epsilon){

		get_pivot(A, N, k, l);

		rotate(A, V, k, l, N);

		iterations += 1;
	}

	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "\nComputation Time = " << time << " s\n";

	cout << "Diagonalized in " << iterations << " iterations\n" << endl;

	cout << "* EIGENVECTORS *" << endl;
	print_matrix(V,N);


	return 0;
}
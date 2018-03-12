#include "jacobi.h"

int main(int argc, char *argv[]){

	int N, k, l;

	if(argc < 1){ cout << "Input number of mesh points N." << endl; }
	else{ 
		N = atoi(argv[1]);
	}	

	double rhomin = 0.0, rhomax = 10.0;
	double h = (rhomax-rhomin)/(N+1), hh = h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double time;

	// set up rho and potential
	vec rho(N), V(N);
	for(int i = 0; i < N; i++){
		rho(i) = rhomin + i*h;
		V(i) = rho(i)*rho(i);
	}

	// set-up Hamiltonian
	mat H = zeros<mat>(N,N);
	for(int i = 0; i < N-1; i++){
		H(i,i) = d+V(i);
		H(i+1,i) = a;
		H(i,i+1) = a;
	}
	H(N-1,N-1) = d+V(N-1);	

	// set-up matrix of eigenvectors
	mat U = eye<mat>(N,N);

	// start timer
	clock_t initial, final;
	initial = clock();

	// solve for eigenvalues and eigenvectors
	jacobi(H, U, N);

	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "Computation Time = " << time << " s\n";

	write_n_eigs(H, U, rho, N, 3, "oscillator");

	return 0;
}
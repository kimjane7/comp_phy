#include "jacobi.h"

int main(int argc, char *argv[]){

	int N, k, l;
	double omega;

	if(argc < 3){ cout << "Input number of mesh points N and at least one freqency omega_r" << endl; }
	else{
		N = atoi(argv[1]);
		for(int iarg = 2; iarg < argc; iarg++){

			omega = atof(argv[iarg]);
		
			double rhomin = 0.0, rhomax = 7.0;
			double h = (rhomax-rhomin)/(N+1), hh = h*h;
			double d = 2.0/hh, a = -1.0/hh;
			double epsilon = 1E-8, time;


			// set up rho and potential, with (V) and without (V0) Coulomb interaction
			vec rho(N), V0(N), V(N);
			for(int i = 0; i < N; i++){
				rho(i) = rhomin + (i+1.0)*h;
				V0(i) = omega*omega*rho(i)*rho(i);
				V(i) = V0(i) + 1.0/rho(i);
			}

			// set-up Hamiltonians, with (H) and without (H0) Coulomb interaction
			mat H0 = zeros<mat>(N,N);
			mat H = zeros<mat>(N,N);
			for(int i = 0; i < N-1; i++){
				H0(i,i) = d+V0(i);
				H0(i+1,i) = a;
				H0(i,i+1) = a;
				H(i,i) = d+V(i);
				H(i+1,i) = a;
				H(i,i+1) = a;
			}
			H(N-1,N-1) = d+V(N-1);
			H0(N-1,N-1) = d+V0(N-1);

			// set-up matrix of eigenvectors
			mat U0 = eye<mat>(N,N);
			mat U = eye<mat>(N,N);

			// start timer
			clock_t initial, final;
			initial = clock();

			// solve for eigenvalues and eigenvectors
			while(offdiag_sq(H0,N) > epsilon){
				get_pivot(H0, N, k, l);
				rotate(H0, U0, k, l, N);
			}

			// solve for eigenvalues and eigenvectors
			while(offdiag_sq(H,N) > epsilon){
				get_pivot(H, N, k, l);
				rotate(H, U, k, l, N);
			}

			// end timer
			final = clock();
			time = (final-initial)/((double) CLOCKS_PER_SEC);
			cout << "\ncomputation time = " << time << " s\n";

			write_ground_state(H, H0, U, U0, rho, N, omega, "dots"+to_string(iarg-1)+".dat");
		}
	} 
	
	return 0;
}
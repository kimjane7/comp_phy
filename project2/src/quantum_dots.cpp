#include "jacobi.h"

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

int main(int argc, char *argv[]){

	int N, k, l;
	double omega;
	string filename;

	if(argc < 4){ cout << "Input number of mesh points N, freqency omega_r, and file name." << endl; }
	else{ 
		N = atoi(argv[1]);
		omega = atof(argv[2]);
		filename = argv[3];
	}	

	double rhomin = 0.0, rhomax = 7.0;
	double h = (rhomax-rhomin)/N, hh = h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double epsilon = 1E-8, time, u, u0;


	// set up rho and potential, with (V) and without (V0) Coulomb interaction
	vec rho(N), V0(N), V(N);
	for(int i = 0; i < N; i++){
		rho(i) = rhomin + (i+0.5)*h;
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

	// write to file
	cout << "writing to file ---> " << filename << endl;
	ofstream ofile;
	ofile.open(filename);
	ofile << "# N = " << N << endl;
	ofile << "# omega_r = " << omega << endl;
	ofile << "# computation time = " << time << endl;
	ofile << "# first eigenvalue [Coulomb] = " << get_nth_eigenvalue(H, N, 0) << endl;
	ofile << "# first eigenvalue [no Coulomb] = " << get_nth_eigenvalue(H0, N, 0) << endl;
	ofile << "# rho, |u(rho)|^2 [Coulomb], |u0(rho)|^2 [no Coulomb], V(rho), V0(rho)" << endl;
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

	cout << "done!\n" << endl;

	return 0;
}
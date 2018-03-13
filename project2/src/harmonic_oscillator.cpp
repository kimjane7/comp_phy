#include "jacobi.h"

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
		ofile << "# N = " << N << endl;
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
		rho(i) = rhomin + (i+1)*h;
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

	write_n_eigs(H, U, rho, N, 6, "oscillator");

	return 0;
}
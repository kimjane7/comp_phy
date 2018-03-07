#include "jacobi.h"

void check_eigenvalues(mat& H, int N, double a, double d){

	vec eigvals(N);
	double lambda;

	for(int i = 0; i < N; i++){
		eigvals(i) = H(i,i);
	}

	sort(eigvals.begin(),eigvals.end());

	// check if eigenvalues are correct
	cout << left << "* EIGENVALUES *\n";
	cout << showpoint;
	cout << setw(15) << "Calculated:";
	cout << setw(15) << "Exact:" << endl;

	for(int i = 0; i < N; i++){

		lambda = 3.0+i*4.0;

		cout << setprecision(10) << setw(15) << eigvals(i); 
		cout << setprecision(10) << setw(15) << lambda << endl;		
	}

	cout << endl;
}

void print_to_file(mat& D, mat& U, vec& rho, int N, string filename){


	ofstream ofile;
	double lambda, u, uu;
	vec eigvals(N);
	ivec index(3);

	// sort eigenvalues by ascending order
	for(int i = 0; i < N; i++){
		eigvals(i) = D(i,i);
	}
	sort(eigvals.begin(),eigvals.end());

	// get indices of lowest 3 eigenvalues
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < N; j++){
			if(eigvals(i) == D(j,j)) index(i) = j;
		}
	}

	// print u(rho) for lowest three eigenvalues
	for(int j = 0; j < 3; j++){

		string outfile = filename + to_string(j+1) + ".dat";
		lambda = 3.0+j*4.0;

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
	}
}


int main(int argc, char *argv[]){

	int N, k, l;

	if(argc < 1){ cout << "Input number of mesh points N." << endl; }
	else{ 
		N = atoi(argv[1]);
	}	

	double rhomin = 0.0, rhomax = 10.0;
	double h = (rhomax-rhomin)/N, hh = h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double epsilon = 1E-8, time;

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
	int iterations = 0;
	while(offdiag_sq(H,N) > epsilon){

		get_pivot(H, N, k, l);

		rotate(H, U, k, l, N);

		iterations += 1;
	}

	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "\nComputation Time = " << time << " s\n";

	cout << "Diagonalized in " << iterations << " iterations\n" << endl;

	//check_eigenvalues(H,N,a,d);

	//cout << "* EIGENVECTORS *" << endl;
	//print_matrix(U,N);

	print_to_file(H, U, rho, N, "oscillator");

	return 0;
}
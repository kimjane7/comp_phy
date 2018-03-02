#include <string> 
#include "jacobi.h"

const double hbar = 1.0545718E-34; // Js

void check_eigenvalues(mat& H, int n, double a, double d){

	vec eigvals(n);
	double lambda;

	for(int i = 0; i < n; i++){
		eigvals(i) = H(i,i);
	}

	sort(eigvals.begin(),eigvals.end());

	// check if eigenvalues are correct
	cout << left << "* EIGENVALUES *\n";
	cout << showpoint;
	cout << setw(15) << "Calculated:";
	cout << setw(15) << "Exact:" << endl;

	for(int i = 0; i < n; i++){

		lambda = 3.0+i*4.0;

		cout << setprecision(10) << setw(15) << eigvals(i); 
		cout << setprecision(10) << setw(15) << lambda << endl;		
	}

	cout << endl;
}

void print(mat& D, mat& U, vec& rho, int n, string filename){


	ofstream ofile;
	double lambda;
	vec eigvals(n);
	ivec index(3);

	// sort eigenvectors by ascending order
	for(int i = 0; i < n; i++){
		eigvals(i) = D(i,i);
	}
	sort(eigvals.begin(),eigvals.end());

	// get indices of lowest 3 eigenvalues
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < n; j++){
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

		for(int i = 0; i < n; i++){
			ofile << rho(i) << "\t" << U(i,index(j)) << endl;
		}
		ofile.close();
	}
}

int main(int argc, char *argv[]){

	int n, k, l;
	double m, w;

	if(argc < 1){ cout << "Input number of mesh points n." << endl; }
	else{ 
		n = atoi(argv[1]);
	}	

	double rhomin = 0.0, rhomax = 5.0;
	double h = (rhomax-rhomin)/n, hh = h*h;
	double d = 2.0/hh, a = -1.0/hh;
	double epsilon = 1E-8, time;

	// set up rho and potential vectors
	vec rho(n), V(n);
	for(int i = 0; i < n; i++){
		rho(i) = rhomin + i*h;
		V(i) = rho(i)*rho(i);
	}

	// set-up Hamiltonian
	mat H = zeros<mat>(n,n);
	for(int i = 0; i < n-1; i++){
		H(i,i) = d+V(i);
		H(i+1,i) = a;
		H(i,i+1) = a;
	}
	H(n-1,n-1) = d+V(n-1);	

	// set-up matrix of eigenvectors
	mat U = eye<mat>(n,n);

	// start timer
	clock_t initial, final;
	initial = clock();

	// solve for eigenvalues and eigenvectors
	int iterations = 0;
	while(offdiag_sq(H,n) > epsilon){

		get_pivot(H, n, k, l);

		rotate(H, U, k, l, n);

		iterations += 1;
	}

	// end timer
	final = clock();
	time = (final-initial)/((double) CLOCKS_PER_SEC);
	cout << "\nComputation Time = " << time << " s\n";

	cout << "Diagonalized in " << iterations << " iterations\n" << endl;

	check_eigenvalues(H,n,a,d);

	//cout << "* EIGENVECTORS *" << endl;
	//print_matrix(U,n);

	print(H, U, rho, n, "oscillator");

	return 0;
}
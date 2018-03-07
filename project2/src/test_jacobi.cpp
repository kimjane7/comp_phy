#include "catch.hpp"
#include "jacobi.h"

TEST_CASE("TEST: Frobenius norm"){

	int N = 3, k, l;

	mat A = zeros<mat>(N,N);
	for(int i = 0; i < N-1; i++){
		A(i,i) = 2.0;
		A(i+1,i) = -1.0;
		A(i,i+1) = -1.0;
	}
	A(N-1,N-1) = 2.0;

	mat U = eye<mat>(N,N);

	double norm_A = norm_sq(A,N);
	double epsilon = 1E-5;

	while(offdiag_sq(A, N) > epsilon){
		get_pivot(A, N, k, l);
		rotate(A, U, k, l, N);
	}

	double norm_D = norm_sq(A,N);

	REQUIRE((norm_A-norm_D)==Approx(0.0));
}

TEST_CASE("TEST: Largest off-diagonal element"){

	int N = 4, k, l;
	mat A = eye<mat>(N,N);
	mat U = eye<mat>(N,N);
	for(int i = 0; i < N-1; i++){
		A(i,i+1) = i+1;
	}

	get_pivot(A, N, k, l);

	REQUIRE(k==N-2);
	REQUIRE(l==N-1);
}

TEST_CASE("TEST: Tridiagonal Toeplitz matrix eigenvalues"){

	int N = 5, k, l;
	double d = 200.0, a = -100.0;
	double epsilon = 1E-5, lambda;
	vec eigvals(N);
	vec exact(N);

	mat A = zeros<mat>(N,N);
	for(int i = 0; i < N-1; i++){
		A(i,i) = d;
		A(i+1,i) = a;
		A(i,i+1) = a;
	}
	A(N-1,N-1) = d;

	mat U = eye<mat>(N,N);
	
	while(offdiag_sq(A, N) > epsilon){
		get_pivot(A, N, k, l);
		rotate(A, U, k, l, N);
	}

	for(int i = 0; i < N; i++){
		eigvals(i) = A(i,i);
		exact(i) = d+2.0*a*cos((i+1)*pi/(N+1));
	}

	sort(eigvals.begin(),eigvals.end());

	for(int i = 0; i < N; i++){
		REQUIRE(eigvals(i)==Approx(exact(i)));
	}
}

TEST_CASE("TEST: Eigenvector orthogonality"){

	int N = 6, k, l;
	double dot, epsilon = 1E-5;;

	mat R = randu<mat>(N,N);
	mat A = R.t()*R;
	mat U = eye<mat>(N,N);

	while(offdiag_sq(A, N) > epsilon){
		get_pivot(A, N, k, l);
		rotate(A, U, k, l, N);
	}

	for(int j1 = 0; j1 < N; j1++){
		for(int j2 = 0; j2 < N; j2++){
			if(j1 != j2){
				dot = 0.0;
				for(int i = 0; i < N; i++){
					dot += U(i,j1)*U(i,j2);
				}
				REQUIRE(dot == Approx(0.0));
			}
		}
	}
}

TEST_CASE("TEST: Compare eigenvalues with Armadillo"){

	int N = 6, k, l;
	double epsilon = 1E-5;

	mat R = randu<mat>(N,N);
	mat A = R.t()*R;
	mat U = eye<mat>(N,N);

	vec eigvals_arma = eig_sym(A);
	vec eigvals(N);
	
	
	while(offdiag_sq(A, N) > epsilon){
		get_pivot(A, N, k, l);
		rotate(A, U, k, l, N);
	}

	for(int i = 0; i < N; i++){
		eigvals(i) = A(i,i);
	}

	sort(eigvals.begin(),eigvals.end());

	for(int i = 0; i < N; i++){
		REQUIRE(eigvals(i) == Approx(eigvals_arma(i)));
	}

}
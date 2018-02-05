project 1:

	benchmark:

		- contains sample file outputs for programs in folder src.
		- each file is named according to corresponding program output and maximum power of 10^i.
		- additional information printed in each file. 

		figures:
			- contains python code for matplotlib figure
			- contains pdf figure 
	report:

		- contains tex and pdf of report

	src:
		-contains code to solve 1D Poisson eqn (-u''(x)=f(x) --> Av=b) with four different assumptions:

			(1) A is dense --> LU_decomp.cpp
			(2) A is tridiagonal --> gen_tridiag.cpp
			(3) A is symmetric, tridiagonal --> symm_tridiag.cpp
			(4) A has 2's along main diagonal, -1's along both secondary diagonals, and 0's elsewhere.

			- each program requires a maximum power of 10 and file name inline
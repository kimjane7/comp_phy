TEST CLONE

project 1:

	report:

		- contains tex and pdf of report
		- contains python code for matplotlib figure
		- contains pdf figure 

		benchmark:

			- contains sample file outputs for programs in folder src.
			- each file is named according to corresponding program output and maximum power of 10^i. e.g. LU3.dat is the output from running LU_decomp.cpp with n = 10^3.
			- additional information printed in each file. 

	src:
		-contains code to solve 1D Poisson eqn (-u''(x)=f(x) --> Av=b) with four different assumptions:

			(1) A is dense --> LU_decomp.cpp
			(2) A is tridiagonal --> gen_tridiag.cpp
			(3) A is symmetric, tridiagonal --> symm_tridiag.cpp
			(4) A has 2's along main diagonal, -1's along both secondary diagonals, and 0's elsewhere.

			- each program requires a maximum power of 10 and file name inline

project 2:

	report:

		- contains tex and pdf of report
		- contains two python codes for figures with corresponding names

		benchmark:

			- files with the format "oscillator(n).dat" (output of harmonic_oscillator.cpp) contains the wavefunctions for nth energy state of one electron in a harmonic oscillator potential.

			- files with the format "dot(n).dat" (output of quantum_dots.cpp) contains the modified ground state wavefunctions for two electrons a harmonic oscillator potential, with and without Coulomb interactions. The index n corresponds to different harmonic oscillator frequencies, which is printed at the top of each file.

	src:

		- jacobi.cpp & jacobi.h contains function declarations and definitions for jacobi's rotation algorithm.

		- toeplitz.cpp implements jacobi's algorithm for buckling beam problem. requires size of matrix N inline.

		- harmonic_oscillator.cpp jacobi's algorithm to solve for wavefunctions and energies of one electron in a harmonic oscillator potential. requires number of mesh points N inline.

		- quantum_dots.cpp implements jacobi's algorithm to solve for wavefunctions and energies of two electrons in a harmonic oscillator potential, with and without Coulomb interactions. requires number of mesh points N AND a list of at least one harmonic oscillator frequency omega_r inline.

		- test_main.cpp & catch.hpp are used for unit testing.

		- test_jacobi.cpp has definitions of unit tests for jacobi's algorithm. 

		- makefile creates executables for toeplitz.cpp, harmonic_oscillator.cpp, quantum_dots.cpp, and test_jacobi.cpp. The executables are names toeplitz, harmonic, quantum, and testcode, respectively. 

		- bisect.cpp & bisect.h contains (unfinished!) function declarations and definitions for a different diagonalization algorithm. 

project 3:

	report:

		- contains tex and pdf of report
		- contains three python codes for plotting figures

		figures:

			- contains orbit and energy plots for each system

		benchmark:

			- file names have the format "system_frame_method{n}_content.dat", where system = {binary, trinary, solar_system}, frame = {fixed, CM}, method = {euler, vv}, n is an integer corresponding to N, and content = {X, Y, E}. 

			- benchmark files for entire solar system was removed because they were too large to push at one point. even after they were made smaller, the files in the history made me unable to push any later changes. 

	src: 

		- main.cpp contains main program

		- makefile contains instructions for compiling program and unit test

		- solar_system.cpp & solar_system.h contain declarations and definitions of CSolarSystem class.

		- planet.cpp & planet.h contain declarations and definitions of CPlanet class.

		- catch.hpp & test_main.cpp are necessary to compile unit tests

		- unit_test.cpp contains definitions of unit tests


project 4:

	report:

		- contains tex and pdf of report
		- contains four python codes for plotting figures

		figures:

			- contains phase portraits and monte carlo samples plots for populations A, B, C, and D.

			- contains results of deterministic approximation

		benchmark:

			- contains S, I, and R for deterministic and monte carlo approaches and all populations

			- numbers correspond to the individual samples

			- if file name contains "stats", it contains the averages and standard deviations of all the samples for that population


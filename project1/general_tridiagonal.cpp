/* Jane
*/


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>


int main(int argc, char *argv[]){

	int max = atoi(argv[1]); // highest power of 10
	int n;					 // dimension of equation
	double h;				 // spacing
	double f;				 // place holder for function value
	double x;				 // place holder for function input


	// loop through powers of 10
	for(int i = min; i <= max; i++){

		n = (int) pow(10.0,i); // is int large enough for n?
		h = 1.0/n;
		hh = h*h;

		double *u = new double[n]; // upper diagonal of A
		double *d = new double[n]; // diagonal of A
		double *l = new double[n]; // lower diagonal of A
		double *b = new double[n]; // Av=b
		double *v = new double[n]; // Av=b
		double *v_exact = new double[n]; // Av=b

		// set-up
		for(int j = 0; j < n; j++){

			// fill out b and v_exact vectors
			x = j*h;
			f = 100.0*exp(-10.0*x);
			b[j] = hh*f;
			v_exact[j] = 1.0-(1.0-exp(-10.0))*x-exp(-10*x);
			    
			// fill out matrix diagonals
			//read from file for general matrix?
			u[j] = -1.0;
			d[j] = 2.0;
			u[j] = -1.0;



		}

	}
	
	



	// momentum (150-2450 MeV in x-direction)
	std::vector<double> p; 
	p.resize(4);
	p[2] = 0.0;
	p[3] = 0.0;
	
	FILE *fptr;
	fptr = fopen(argv[2], "w");
	fprintf(fptr,"# p_x (MeV)    E*dN/d^3p\n");

	// sweep through different values of p_x
	for(p[1] = 110.0 ; p[1] < 2500.0; p[1] += 10.0) { 

		// energy
		p[0] = pow(p[1]*p[1]+m*m,0.5);

		//integral with respect to tau
		EdN_d3p = 0.0;
		for(tau = tau0+0.5*dtau; tau < tauf; tau += dtau) {

			// radius of surface
			R = R0-(tau-tau0)*v;

			// radial velocity
			u = R/R0;

			// bessel functions used in eta integral
			K0 = boost::math::cyl_bessel_k(0, pow((1.0+u*u),0.5)*p[0]/T );
			K1 = boost::math::cyl_bessel_k(1, pow((1.0+u*u),0.5)*p[0]/T );
			I0 = boost::math::cyl_bessel_i(0, p[1]*u/T );
			I1 = boost::math::cyl_bessel_i(1, p[1]*u/T );


			EdN_d3p += (p[0]*v*K1*I0+p[1]*K0*I1)*R*tau*dtau;
		}

		EdN_d3p = 4*pi*prefactor*EdN_d3p;

		fprintf(fptr,"%lf\t%.10e\n",p[1],EdN_d3p);
		
	}

	fclose(fptr);

	return 0;
}
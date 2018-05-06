#include "SIRS.h"

int main(int argc, char *argv[]){

	CInfectedPopulation A(400, 4.0, 1.0, 0.5);
	A.lattice_SIRS("lattice_A", 100, 300, 100, 20.0);
	/*
	A.montecarlo_SIRS("montecarlo_A", 100, 300, 100, 20.0);
	A.deterministic_SIRS("deterministic_A", 300.0, 100, 20.0);
	A.generate_phaseportrait("phaseportrait_A", 15.0);
	
	CInfectedPopulation B(400, 4.0, 2.0, 0.5);
	B.montecarlo_SIRS("montecarlo_B", 100, 300, 100, 20.0);
	B.deterministic_SIRS("deterministic_B", 300.0, 100, 20.0);
	B.generate_phaseportrait("phaseportrait_B", 15.0);


	CInfectedPopulation C(400, 4.0, 3.0, 0.5);
	C.montecarlo_SIRS("montecarlo_C", 100, 300, 100, 20.0);
	C.deterministic_SIRS("deterministic_C", 300.0, 100, 20.0);
	C.generate_phaseportrait("phaseportrait_C", 15.0);

	CInfectedPopulation D(400, 4.0, 4.0, 0.5);
	D.montecarlo_SIRS("montecarlo_D", 100, 300, 100, 20.0);
	D.deterministic_SIRS("deterministic_D", 300.0, 100, 20.0);
	D.generate_phaseportrait("phaseportrait_D", 15.0);
	*/

	return 0;
}
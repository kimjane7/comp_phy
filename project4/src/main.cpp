#include "SIRS.h"

int main(int argc, char *argv[]){

	CInfectedPopulation A(500, 5.0, 0.5, 2.0);
	A.montecarlo_SIRS("montecarlo_A", 100, 400, 100, 15.0);
	A.deterministic_SIRS("deterministic_A", 400.0, 100, 15.0);
	A.generate_phaseportrait("phaseportrait_A", 15.0);

	CInfectedPopulation B(500, 1.5, 0.5, 2.0);
	B.montecarlo_SIRS("montecarlo_B", 100, 400, 100, 15.0);
	B.deterministic_SIRS("deterministic_B", 400.0, 100, 15.0);
	B.generate_phaseportrait("phaseportrait_B", 15.0);

	return 0;
}
#include "SIRS.h"

int main(int argc, char *argv[]){

	CInfectedPopulation test(1000, 5.0, 1.0, 2.0);

	test.deterministic_SIRS("deterministic_test.dat", 900.0, 100.0, 0.0, 10.0);
	test.stochastic_SIRS("stochastic_test.dat", 900, 100, 0, 10.0);

	return 0;
}
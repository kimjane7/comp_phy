#include "SIRS.h"

int main(int argc, char *argv[]){

	CInfectedPopulation test(500, 5.0, 0.5, 2.0);

	test.deterministic_SIRS("deterministic_test", 400.0, 100.0, 15.0);
	test.stochastic_SIRS("stochastic_test", 100, 400, 100, 15.0);

	return 0;
}
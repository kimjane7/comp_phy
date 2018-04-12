#include "SIRS.h"

int main(int argc, char *argv[]){

	CDiseasedPopulation test(1000, 5.0, 1.0, 2.0);

	test.deterministic_SIRS("test.dat", 900, 100, 0, 10.0);

	return 0;
}
#ifndef SIRS_H
#define SIRS_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <vector>
#include <random>

using namespace std;

class CInfectedPopulation{
public:

	int N_;
	double dt_;
	double a_, b_, c_;

	CInfectedPopulation();
	CInfectedPopulation(int N, double a, double b, double c);
	~CInfectedPopulation(){};

	void deterministic_SIRS(string filename, double S0, double I0, double tf);
	void generate_phaseportrait(string filename, double tf);
	void montecarlo_SIRS(string filename, int nsamples, int S0, int I0, double tf);
	void lattice_SIRS(string filename, int nsamples, int SO, int I0, double tf);
};

#endif
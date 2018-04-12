#ifndef SIRS_H
#define SIRS_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <boost/random.hpp>

using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::setprecision;
using std::setw;
using std::left;

using boost::uniform_01;
using boost::mt19937;

double random01(mt19937 generator){
	static uniform_01<mt19937> dist(generator);
	return dist();
}

// for now N = constant
class CInfectedPopulation{
public:

	int N_;

	double dt_;
	double a_, b_, c_;

	CInfectedPopulation();
	CInfectedPopulation(int N, double a, double b, double c);
	~CInfectedPopulation(){};

	void deterministic_SIRS(string filename, double S0, double I0, double R0, double tf);
	void stochastic_SIRS(string filename, int S0, int I0, int R0, double tf);
};

#endif
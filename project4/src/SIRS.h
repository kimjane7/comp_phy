#ifndef SIRS_H
#define SIRS_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <vector>
#include "time.h"

using namespace std;


// for now N = constant
class CDiseasedPopulation{
public:

	int N_;

	double dt_;
	double a_, b_, c_;

	CDiseasedPopulation();
	CDiseasedPopulation(int N, double a, double b, double c);
	~CDiseasedPopulation(){};

	void deterministic_SIRS(string filename, double S0, double I0, double R0, double tf);
	void stochastic_SIRS(string filename, int S0, int I0, int R0, double tf);
};

#endif
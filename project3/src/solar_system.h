// solves system of the form:      
//      x'' = fx(x,y)
//      y'' = fy(x,y)

#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <vector>
#include "time.h"
#include "planet.h"

using namespace std;

class CSolarSystem{

public:

	int N_, planets_;
	double t0_, tf_, h_;
	vector<CPlanet> planet_list_;

	CSolarSystem();
	CSolarSystem(int N, double t0, double tf);
	~CSolarSystem(){}

	double fx(double x, double y){ return -4.0*pi*pi*x/pow(x*x+y*y,1.5); }
	double fy(double x, double y){ return -4.0*pi*pi*y/pow(x*x+y*y,1.5); }

	void add(CPlanet NewPlanet);

	void solve_euler(string filename);
	void compare_euler(string filename, int maxpower);

};

#endif
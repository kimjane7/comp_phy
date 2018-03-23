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
	double h_;
	vector<CPlanet> planet_list_;

	vector<double> t_;
	vector<vector<double>> x_, y_, vx_, vy_;

	CSolarSystem();
	CSolarSystem(int N, double t0, double tf);
	~CSolarSystem(){}

	void add(CPlanet NewPlanet);

	double distance(int i, int j, int k);

	double fx(double x, double y){ return -4.0*pi*pi*x/pow(x*x+y*y,1.5); }
	double fy(double x, double y){ return -4.0*pi*pi*y/pow(x*x+y*y,1.5); }

	double ax(int i, int j);
	double ay(int i, int j);

	void solve_euler(string filename);
	void compare_euler(string filename, int maxpower);

	void solve_vv(string filename);
	void compare_vv(string filename, int maxpower);

};

#endif
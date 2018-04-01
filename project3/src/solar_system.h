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

const double pi = 4.0*atan(1.0);
const double prefactor = 4.0*pi*pi;

class CSolarSystem{

public:

	int N_, planets_;
	double h_;
	vector<CPlanet> planet_list_;

	double xCM_, yCM_;
	bool CM_frame_;

	vector<double> t_;
	vector<vector<double>> x_, y_, vx_, vy_;

	CSolarSystem();
	CSolarSystem(int N, double t0, double tf, bool CM_frame);
	~CSolarSystem(){}

	void add(CPlanet NewPlanet);
	void change_N(int new_N);

	double distance(int i, int j, int k);

	void get_acceleration(int i, int j, double& ax, double& ay);
	void get_energy(int i, int j, double& KE, double& PE);

	void initialize();
	void solve_euler();
	void solve_vv();

	void compare_euler(string systemname, int maxpower);
	void compare_vv(string systemname, int maxpower);

	void write_orbits(string filename);
	void write_energies(string filename);
};

#endif
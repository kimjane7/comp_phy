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
#include "time.h"

using namespace std;

const double pi = 4.0*atan(1.0);
const double prefactor = 4.0*pi*pi;

class CSolarSystem{

public:

	bool CM_frame_;
	int N_, planets_, dim_;
	double h_, comp_time_;
	double CM_[3];


	vector<CPlanet> planet_list_;
	vector<double> t_;
	vector<vector<double>> x_, y_, z_;
	vector<vector<double>> vx_, vy_, vz_;

	CSolarSystem();
	CSolarSystem(int N, int dim, double t0, double tf, bool CM_frame);
	~CSolarSystem(){}

	void add(CPlanet NewPlanet);
	void change_N(int new_N);

	double distance(int i, int j, int k);
	double velocity(int i, int j);

	void get_acceleration(int i, int j, double& ax, double& ay, double& az);
	void get_energy(int i, int j, double& E);
	void get_angmomentum(int i, int j, double& L);

	void initialize();

	void solve_euler();
	void solve_vv();

	void write_orbits(string filename);
	void write_energies(string filename);

	void compare_euler(string systemname, int maxpower);
	void compare_vv(string systemname, int maxpower);
};

#endif
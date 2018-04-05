#ifndef PLANET_H
#define PLANET_H

#include <string> 
#include <iomanip>
#include <vector>

using namespace std;

class CPlanet{
public:

	string name_;
	double m_;
	double x0_[3], v0_[3];

	CPlanet();
	CPlanet(string name, double m_ratio, double x0, double y0, double z0, double vx0, double vy0, double vz0);

	void change_IC(double x0, double y0, double z0, double vx0, double vy0, double vz0);

};

#endif
#ifndef PLANET_H
#define PLANET_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <vector>
#include "time.h"


using namespace std;

const double pi = 4.0*atan(1.0);

class CPlanet{
public:

	double m_;
	double x0_[3]. x_[3];
	double v0_[3], v_[3];

	CPlanet();
	CPlanet(double m_ratio, double x0, double y0, double z0, double vx0, double vy0, double vz0);

};

#endif
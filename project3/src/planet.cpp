#include "planet.h"

CPlanet::CPlanet(){

	name_ = "PlanetX"
	m_ = 1.0;

	x0_[0] = 1.0;
	x0_[1] = 0.0;
	x0_[2] = 0.0;

	v0_[0] = 0.0;
	v0_[1] = 0.0;
	v0_[2] = 0.0;

}	

CPlanet::CPlanet(string name, double m_ratio, double x0, double y0, double z0, double vx0, double vy0, double vz0){

	name_ = name;
	m_ = m_ratio;

	x0_[0] = x0;
	x0_[1] = y0;
	x0_[2] = z0;

	v0_[0] = vx0;
	v0_[1] = vy0;
	v0_[2] = vz0;

}

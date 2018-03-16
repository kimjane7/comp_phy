#include "planet.h"

CPlanet::CPlanet(){

	m_ = 1.0;

	x_[0] = 1.0;
	x_[1] = 0.0;
	x_[2] = 0.0;

	v_[0] = 0.0;
	v_[1] = 0.0;
	v_[2] = 0.0;

}	

CPlanet::CPlanet(double m_ratio, double x0, double y0, double z0, double vx0, double vy0, double vz0){

	m_ = m_ratio;

	x_[0] = x0;
	x_[1] = y0;
	x_[2] = z0;

	v_[0] = vx0;
	v_[1] = vy0;
	v_[2] = vz0;

}

double CPlanet::distance(CPlanet OtherPlanet){

	double x, y, z;

	x = (this->x_[0])-(OtherPlanet.x_[0]);
	y = (this->x_[1])-(OtherPlanet.x_[1]);
	z = (this->x_[2])-(OtherPlanet.x_[2]);

	return sqrt(x*x+y*y+z*z);

}

double CPlanet::force(CPlanet OtherPlanet){

	double r = this->distance(OtherPlanet);

	if(r != 0){
		double m1 = this->m_;
		double m2 = OtherPlanet.m_;
		return 4.0*pi*pi*m1*m2/(r*r);
	}

	else return 0.0;
}

#include "planet.h"
#include "solar_system.h"

int main(int argc, char *argv[]){

	// CPlanet planet(name, mass ratio, x0, y0, z0, vx0, vy0, vz0);
	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, 2*pi, 0.0);
	CPlanet jupiter("Jupiter", 9.5E-4, 5.2, 0.0, 0.0, 0.0, 2*pi/sqrt(5.2), 0.0);
	CPlanet mars("Mars", 3.3E-7, 1.52, 0.0, 0.0, 0.0, 2*pi/sqrt(1.52), 0.0);
	CPlanet venus("Venus", 2.45E-6, 0.72, 0.0, 0.0, 0.0, 2*pi/sqrt(0.72), 0.0);
	CPlanet saturn("Saturn", 0.00095, 5.2, 0.0, 0.0, 0.0, 2*pi/sqrt(5.2), 0.0);
	CPlanet mercury("Mercury", 0.00095, 5.2, 0.0, 0.0, 0.0, 2*pi/sqrt(5.2), 0.0);
	CPlanet uranus("Uranus", 0.00095, 5.2, 0.0, 0.0, 0.0, 2*pi/sqrt(5.2), 0.0);
	CPlanet neptune("Neptune", 0.00095, 5.2, 0.0, 0.0, 0.0, 2*pi/sqrt(5.2), 0.0);
	CPlanet pluto("Pluto", 0.00095, 5.2, 0.0, 0.0, 0.0, 2*pi/sqrt(5.2), 0.0);


	// CSolarSystem system(N, t0, tf);
	CSolarSystem trinary_system(1, 0.0, 15.0);
	trinary_system.add(sun);
	trinary_system.add(earth);
	trinary_system.add(jupiter);

	trinary_system.compare_euler("trinary", 6);
	trinary_system.compare_vv("trinary", 4);

	/*
	CSolarSystem binary_system(1, 0.0, 10.0);
	binary_system.add(sun);
	binary_system.add(earth);

	binary_system.compare_euler("binary", 6);
	binary_system.compare_vv("binary", 4);
	*/



	return 0;
}

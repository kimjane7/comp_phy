#include "planet.h"
#include "solar_system.h"

int main(int argc, char *argv[]){

	// CPlanet planet(name, mass ratio, x0, y0, z0, vx0, vy0, vz0);
	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, 2*pi, 0.0);
	CPlanet jupiter("Jupiter", 9.5E-4, 5.2, 0.0, 0.0, 0.0, 2*pi/sqrt(5.2), 0.0);
	CPlanet mars("Mars", 3.3E-7, 1.52, 0.0, 0.0, 0.0, 2*pi/sqrt(1.52), 0.0);
	CPlanet venus("Venus", 2.45E-6, 0.72, 0.0, 0.0, 0.0, 2*pi/sqrt(0.72), 0.0);
	CPlanet saturn("Saturn", 2.75E-4, 9.54, 0.0, 0.0, 0.0, 2*pi/sqrt(9.54), 0.0);
	CPlanet mercury("Mercury", 1.65E-7, 0.39, 0.0, 0.0, 0.0, 2*pi/sqrt(0.39), 0.0);
	CPlanet uranus("Uranus", 4.4E-5, 19.19, 0.0, 0.0, 0.0, 2*pi/sqrt(19.19), 0.0);
	CPlanet neptune("Neptune", 5.15E-5, 30.06, 0.0, 0.0, 0.0, 2*pi/sqrt(30.06), 0.0);
	CPlanet pluto("Pluto", 6.55E-9, 39.53, 0.0, 0.0, 0.0, 2*pi/sqrt(39.53), 0.0);


	// CSolarSystem system(N, t0, tf, CM frame on);
	CSolarSystem solar_system(2E4, 0.0, 1.0, true);
	solar_system.add(sun);
	solar_system.add(earth);
	solar_system.add(jupiter);
	solar_system.add(mars);
	solar_system.add(venus);
	solar_system.add(saturn);
	solar_system.add(mercury);
	solar_system.add(uranus);
	solar_system.add(neptune);
	solar_system.add(pluto);

	solar_system.initialize();
	solar_system.solve_euler();
	solar_system.write_orbits("CM_euler_");
	solar_system.solve_vv();
	solar_system.write_orbits("CM_vv_");

	return 0;
}

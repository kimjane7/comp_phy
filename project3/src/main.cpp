#include "planet.h"
#include "solar_system.h"

int main(int argc, char *argv[]){

	// CPlanet planet(mass ratio, x0, y0, z0, vx0, vy0, vz0);
	CPlanet sun(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth(0.000003, 1.0, 0.0, 0.0, 0.0, 2*pi, 0.0);

	// CSolarSystem system(N, t0, tf);
	CSolarSystem binary_system(1E6, 0.0, 10.0);
	binary_system.add(sun);
	binary_system.add(earth);

	binary_system.compare_euler("binary", 6);

	return 0;
}

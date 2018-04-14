#include "planet.h"
#include "solar_system.h"

int main(int argc, char *argv[]){

	
	// CPlanet planet(name, mass ratio, x0, y0, z0, vx0, vy0, vz0);
	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet mercury("Mercury", 1.65E-7, 0.39, 0.0, 0.0, 0.0, 2.0*pi/sqrt(0.39), 0.0);
	CPlanet venus("Venus", 2.45E-6, 0.72, 0.0, 0.0, 0.0, 2.0*pi/sqrt(0.72), 0.0);
	CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, 2.0*pi, 0.0);
	CPlanet mars("Mars", 3.3E-7, 1.52, 0.0, 0.0, 0.0, 2.0*pi/sqrt(1.52), 0.0);
	CPlanet jupiter("Jupiter", 9.5E-4, 5.2, 0.0, 0.0, 0.0, 2.0*pi/sqrt(5.2), 0.0);
	CPlanet saturn("Saturn", 2.75E-4, 9.54, 0.0, 0.0, 0.0, 2.0*pi/sqrt(9.54), 0.0);
	CPlanet uranus("Uranus", 4.4E-5, 19.19, 0.0, 0.0, 0.0, 2.0*pi/sqrt(19.19), 0.0);
	CPlanet neptune("Neptune", 5.15E-5, 30.06, 0.0, 0.0, 0.0, 2.0*pi/sqrt(30.06), 0.0);
	CPlanet pluto("Pluto", 6.55E-9, 39.53, 0.0, 0.0, 0.0, 2.0*pi/sqrt(39.53), 0.0);


	// CSolarSystem system(N, dim, t0, tf, CM frame on);
	CSolarSystem binary(1, 2, 0.0, 50.0, false);
	binary.add(sun);
	binary.add(earth);
	binary.compare_euler("binary_fixed",7);
	binary.compare_vv("binary_fixed",7);
	

	CSolarSystem trinary(1, 2, 0.0, 60.0, false);
	trinary.add(sun);
	trinary.add(earth);
	trinary.add(jupiter);
	trinary.compare_euler("trinary_fixed",7);
	trinary.compare_vv("trinary_fixed",7);

	/*
	CSolarSystem solar_system(3E6, 2, 0.0, 248.0, false);
	solar_system.add(sun);
	solar_system.add(mercury);
	solar_system.add(venus);
	solar_system.add(earth);
	solar_system.add(mars);
	solar_system.add(jupiter);
	solar_system.add(saturn);
	solar_system.add(uranus);
	solar_system.add(neptune);
	solar_system.add(pluto);
	solar_system.initialize();
	solar_system.solve_euler();
	solar_system.write_orbits("solarsystem_fixed_euler");
	solar_system.solve_vv();
	solar_system.write_orbits("solarsystem_fixed_vv");
	*/
	/*
	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet mercury("Mercury", 1.65E-7, -0.30627, -0.22875, 1.4187E-2, 3.4427, -8.1902, -0.98537);
	CPlanet venus("Venus", 2.45E-6, -0.61163, 0.37715, 4.0389E-2, -3.8569, -6.3546, 0.13531);
	CPlanet earth("Earth", 3.0E-6, 0.93599, 0.35658, -1.4878E-4, -2.3110, 5.8585, -4.2113E-4);
	CPlanet mars("Mars", 3.3E-7, -1.5501, 0.60823, 5.0592E-2, -1.6597, -4.3261, -4.9968E-2);
	CPlanet jupiter("Jupiter", 9.5E-4, -4.5966, -2.9031, 0.11485, 1.4393, -2.1990, 2.3052E-2);
	CPlanet saturn("Saturn", 2.75E-4, -0.36854, -10.049, 0.18938, 1.9241, -0.08120, -0.07517);
	CPlanet uranus("Uranus", 4.4E-5, 17.865, 8.79979, -0.19876, -0.64527, 1.2217, 0.012886);
	CPlanet neptune("Neptune", 5.15E-5, 28.611, -8.8304, -0.47753, 0.33041, 1.1024, -0.030149);
	CPlanet pluto("Pluto", 6.55E-9, 10.538, -31.714, 0.34537, 1.1162, 0.12426, -0.33748);

	CSolarSystem solar_system(3E6, 3, 0.0, 248,  true);
	solar_system.add(sun);
	solar_system.add(mercury);
	solar_system.add(venus);
	solar_system.add(earth);
	solar_system.add(mars);
	solar_system.add(jupiter);
	solar_system.add(saturn);
	solar_system.add(uranus);
	solar_system.add(neptune);
	solar_system.add(pluto);

	solar_system.initialize();
	solar_system.solve_euler();
	solar_system.write_orbits("solarsystem_CM_euler");
	solar_system.solve_vv();
	solar_system.write_orbits("solarsystem_CM_vv");
	*/

	return 0;
}

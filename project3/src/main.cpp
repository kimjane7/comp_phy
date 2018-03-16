#include "planet.h"
#include "solar_system.h"

int main(int argc, char *argv[]){

	CPlanet sun(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth(0.000003, 1.0, 0.0, 0.0, 0.0, 2.0*pi, 0.0);

	CSolarSystem binary_system;
	binary_system.add(sun);
	binary_system.add(earth);

	for(int i = 0; i < binary_system.planets_; i++){
		cout << binary_system.planet_list_[i].m_ << endl;
	}


	return 0;
}
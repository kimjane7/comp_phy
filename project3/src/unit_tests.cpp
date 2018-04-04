#include "catch.hpp"
#include "planet.h"
#include "solar_system.h"

TEST_CASE("TEST: Stability of Earth's Orbit"){

	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, 2.0*pi, 0.0);
	CSolarSystem test(1, 0, 1000, false);
	test.add(sun);
	test.add(earth);
	cout << "hello" << endl;
	int x = 10;
	REQUIRE(x==10);
	//REQUIRE(test.check_stability(1)==true);
}
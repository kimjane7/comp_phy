#include "catch.hpp"
#include "planet.h"
#include "solar_system.h"


TEST_CASE("TEST: Stability of Earth's Orbit"){

	cout << "TEST: Stability of Earth's Orbit" << endl;

	int N = 1E7;

	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, 2.0*pi, 0.0);
	CSolarSystem test(N, 2, 0, 1000, false);
	test.add(sun);
	test.add(earth);
	test.initialize();
	test.solve_vv();

	// check distance between Sun and Earth is 1 AU throughout 1000 years
	for(int i = 0; i <= N; i += 1000) REQUIRE(test.distance(i,0,1) == Approx(1.0));
}

TEST_CASE("TEST: Conservation of Energy"){

	cout << "TEST: Conservation of Energy" << endl;

	int N = 1E4;
	double E0, E;

	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, 2.0*pi, 0.0);
	CSolarSystem test(N, 2, 0, 1, false);
	test.add(sun);
	test.add(earth);
	test.initialize();
	test.solve_vv();

	// get initial energy of Earth
	test.get_energy(0, 1, E0);

	// check energy of Earth remains the same
	for(int i = 0; i <= N; i += 100){
		test.get_energy(i, 1, E);
		REQUIRE(E == Approx(E0));
	}
}

TEST_CASE("TEST: Conservation of Angular Momentum"){

	cout << "TEST: Conservation of Angular Momentum" << endl;

	int N = 1E4;
	double L0, L;

	CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, 2.0*pi, 0.0);
	CSolarSystem test(N, 2, 0, 1, false);
	test.add(sun);
	test.add(earth);
	test.initialize();
	test.solve_vv();

	// get initial energy of Earth
	test.get_angmomentum(0, 1, L0);

	// check energy of Earth remains the same
	for(int i = 0; i <= N; i += 100){
		test.get_angmomentum(i, 1, L);
		REQUIRE(L == Approx(L0));
	}
}

TEST_CASE("TEST: Escape Velocity of Earth"){

	cout << "TEST: Escape Velocity" << endl;

	// bisection method to find escape velocity
	int N = 1E5;
	double a = 2.0*pi, b = 3.0*pi;
	double v0, E;
	while(b-a > 1.0E-5){

		v0 = 0.5*(a+b);

		CPlanet sun("Sun", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		CPlanet earth("Earth", 3.0E-6, 1.0, 0.0, 0.0, 0.0, v0, 0.0);
		CSolarSystem test(N, 2, 0, 1, false);
		test.add(sun);
		test.add(earth);
		test.initialize();
		test.solve_vv();
		test.get_energy(N, 1, E);

		if(E < 0) a = v0;
		if(E > 0) b = v0;
	}

	cout << "\tEstimated escape velocity = " << v0 << " AU/yr." << endl;
	cout << "\tTheoretical escape velocity = " << 2.0*sqrt(2.0)*pi << " AU/yr." << endl;

	REQUIRE(v0 == Approx(2.0*sqrt(2.0)*pi));
}
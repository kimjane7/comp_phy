#include "solar_system.h"

CSolarSystem::CSolarSystem(){

	N_ = 1000;
	planets_ = 0;
	h_ = 10.0/N_;

	// set-up time vector
	t_.resize(N_+1);
	t_[0] = 0.0;
	t_[N_] = 10.0;
	for(int i = 1; i < N_; i++){ t_[i] = i*h_; }

	// set-up matrices to store data
	x_.resize(N_+1);
	y_.resize(N_+1);
	vx_.resize(N_+1);
	vy_.resize(N_+1);
}

CSolarSystem::CSolarSystem(int N, double t0, double tf){

	N_ = N;
	planets_ = 0;
	h_ = (tf-t0)/N_;

	// set-up time vector
	t_.resize(N_+1);
	t_[0] = t0;
	t_[N_] = tf;
	for(int i = 1; i < N_; i++){ t_[i] = t0+i*h_; }

	// set-up matrices to store data
	x_.resize(N_+1);
	y_.resize(N_+1);
	vx_.resize(N_+1);
	vy_.resize(N_+1);
}

void CSolarSystem::add(CPlanet NewPlanet){

	planets_ += 1;
	planet_list_.push_back(NewPlanet);

	// add one column for each body in the solar system
	for(int i = 0; i < N_+1; i++){
		x_[i].resize(planets_);
		y_[i].resize(planets_);
		vx_[i].resize(planets_);
		vy_[i].resize(planets_);
	}
}

// changes the number of mesh points N and everything dependent on N
// t0, tf, and number of planets in system remain fixed
void CSolarSystem::change_N(int new_N){

	// keep initial and final times
	double t0 = t_[0];
	double tf = t_[N_];

	N_ = new_N;
	h_ = (tf-t0)/N_;

	// update time vector
	t_.resize(N_+1);
	t_[0] = t0;
	t_[N_] = tf;
	for(int i = 1; i < N_; i++){ t_[i] = t0+i*h_; }

	// update matrices to store data
	x_.resize(N_+1);
	y_.resize(N_+1);
	vx_.resize(N_+1);
	vy_.resize(N_+1);

	// update number of columns
	for(int i = 0; i <= N_; i++){
		x_[i].resize(planets_);
		y_[i].resize(planets_);
		vx_[i].resize(planets_);
		vy_[i].resize(planets_);
	}
}

// calculates distance between jth and kth planets at ith time step
double CSolarSystem::distance(int i, int j, int k){

	double x, y;

	x = x_[i][j]-x_[i][k];
	y = y_[i][j]-y_[i][k];

	return sqrt(x*x+y*y);
}

// calculates the acceleration of the jth planet due to 
// the other planets in the system at the ith time step
void CSolarSystem::get_acceleration(int i, int j, double& ax, double& ay){

	double m, r, rrr;

	ax = 0.0;
	ay = 0.0;
	for(int k = 0; k < planets_; k++){
		if(k != j) {
			m = planet_list_[k].m_;
			r = distance(i, j, k);
			rrr = r*r*r;
			// pretty sure these have the right sign, check later
			ax += prefactor*m*(x_[i][k]-x_[i][j])/rrr;
			ay += prefactor*m*(y_[i][k]-y_[i][j])/rrr;		
		}
	}	
}

// calculates the energy of the jth planet at the ith time step
void CSolarSystem::get_energy(int i, int j, double& KE, double& PE){

	double M, m, r;

	M = planet_list_[j].m_;
	KE = 0.5*M*(vx_[i][j]*vx_[i][j]+vy_[i][j]*vy_[i][j]);
	PE = 0.0;
	for(int k = 0; k < planets_; k++){
		if(k != j){
			m = planet_list_[k].m_;
			r = distance(i, j, k);
			PE += -prefactor*M*m/r;
		}
	}
}

// solves for trajectories and energies of all planets
// for now, sun is FIXED at (0,0)
void CSolarSystem::solve_euler(){

	double m, ax, ay;

	// initial conditions of planets (including sun)
	for(int j = 0; j < planets_; j++){
		x_[0][j] = planet_list_[j].x0_[0];
		y_[0][j] = planet_list_[j].x0_[1];
		vx_[0][j] = planet_list_[j].v0_[0];
		vy_[0][j] = planet_list_[j].v0_[1];
	}

	// forward euler
	for(int i = 1; i <= N_; i++){

		// fix the sun at (0,0)
		x_[i][0] = 0.0;
		y_[i][0] = 0.0;
		vx_[i][0] = 0.0;
		vy_[i][0] = 0.0;

		// calculate orbits of other planets
		for(int j = 1; j < planets_; j++){

			x_[i][j] = x_[i-1][j] + h_*vx_[i-1][j];
			y_[i][j] = y_[i-1][j] + h_*vy_[i-1][j];
			get_acceleration(i-1, j, ax, ay);
			vx_[i][j] = vx_[i-1][j] + h_*ax;
			vy_[i][j] = vy_[i-1][j] + h_*ay;
		}
	}
}

void CSolarSystem::solve_vv(){

	double m, ax1, ay1, ax2, ay2;

	// initial conditions of planets (including sun)
	for(int j = 0; j < planets_; j++){
		x_[0][j] = planet_list_[j].x0_[0];
		y_[0][j] = planet_list_[j].x0_[1];
		vx_[0][j] = planet_list_[j].v0_[0];
		vy_[0][j] = planet_list_[j].v0_[1];
	}

	// velocity verlet
	for(int i = 1; i <= N_; i++){

		// fix the sun at (0,0)
		x_[i][0] = 0.0;
		y_[i][0] = 0.0;
		vx_[i][0] = 0.0;
		vy_[i][0] = 0.0;

		// calculate orbits of other planets
		for(int j = 1; j < planets_; j++){

			get_acceleration(i-1, j, ax1, ay1);
			x_[i][j] = x_[i-1][j] + h_*vx_[i-1][j] + 0.5*h_*h_*ax1;
			y_[i][j] = y_[i-1][j] + h_*vy_[i-1][j] + 0.5*h_*h_*ay1;

			get_acceleration(i, j, ax2, ay2);
			vx_[i][j] = vx_[i-1][j] + 0.5*h_*(ax1+ax2);
			vy_[i][j] = vy_[i-1][j] + 0.5*h_*(ay1+ay2);
		}		
	}
}

void CSolarSystem::write_orbits(string filename){

	string Xfile = filename + "X.dat";
	string Yfile = filename + "Y.dat";
	
	// make heading to label columns of file
	string heading = "# time";
	for(int j = 0; j < planets_; j++){
		heading += ", " + planet_list_[j].name_;
	}

	cout << "Writing x-coordinates to '" << Xfile << "'... " << endl;
	cout << "Writing y-coordinates to '" << Yfile << "'... " << endl;

	ofstream outX, outY;
	outX.open(Xfile);
	outY.open(Yfile);

	outX << "# N = " << N_ << endl;
	outX << "# h = " << h_ << endl;
	outX << "# t0 = " << t_[0] << endl;
	outX << "# tf = " << t_[N_] << endl;

	outY << "# N = " << N_ << endl;
	outY << "# h = " << h_ << endl;
	outY << "# t0 = " << t_[0] << endl;
	outY << "# tf = " << t_[N_] << endl;

	outX << heading << endl;
	outY << heading << endl;

	// write positions of all planets to file
	for(int i = 0; i <= N_; i++){

		// time
		outX << left << setw(14) << setprecision(7) << t_[i];
		outY << left << setw(14) << setprecision(7) << t_[i];

		// positions of planets
		for(int j = 0; j < planets_; j++){
			outX << left << setw(14) << setprecision(7) << x_[i][j];
			outY << left << setw(14) << setprecision(7) << y_[i][j];
		}

		outX << endl;
		outY << endl;
	}

	outX.close();
	outY.close();
}

void CSolarSystem::write_energies(string filename){

	double KE, PE;
	string Efile = filename + "E.dat";

	// make heading to label columns of file
	string heading = "# time";
	for(int j = 0; j < planets_; j++){
		heading += ", " + planet_list_[j].name_;
	}

	cout << "Writing energies to '" << Efile << "'... " << endl;

	ofstream outE;
	outE.open(Efile);

	outE << "# N = " << N_ << endl;
	outE << "# h = " << h_ << endl;
	outE << "# t0 = " << t_[0] << endl;
	outE << "# tf = " << t_[N_] << endl;

	outE << heading << endl;

	// write all energies to file
	for(int i = 0; i <= N_; i++){

		// time
		outE << left << setw(14) << setprecision(7) << t_[i];

		// energies of planet
		for(int j = 0; j < planets_; j++){
			get_energy(i, j, KE, PE);
			outE << left << setw(14) << setprecision(7) << KE+PE;
		}

		outE << endl;
	}

	outE.close();
}

// writes orbits and energies to files for different N using euler's method
void CSolarSystem::compare_euler(string systemname, int maxpower){

	cout << "\n*** FORWARD EULER ***\n" << endl;
	for(int n = 2; n <= maxpower; n++){

		cout << "N = 10E" << n << endl;
		string filename = systemname + "_euler_" + to_string(n);
		change_N(pow(10,n));
		solve_euler();
		write_orbits(filename);
		write_energies(filename);
		cout << endl;		
	}
}

// writes orbits and energies to files for different N using velocity verlet method
void CSolarSystem::compare_vv(string systemname, int maxpower){

	cout << "\n*** VELOCITY VERLET ***\n" << endl;
	for(int n = 2; n <= maxpower; n++){

		cout << "N = 10E" << n << endl;
		string filename = systemname + "_vv_" + to_string(n);
		change_N(pow(10,n));
		solve_vv();
		write_orbits(filename);
		write_energies(filename);
		cout << endl;		
	}
}

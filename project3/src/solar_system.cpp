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
	t_.resize(N_+1)
	t_[0] = t0;
	t_[N_] = tf;
	for(int i = 1; i < N_; i++){ t[i] = t0+i*h_; }

	// set-up matrices to store data
	x_.resize(N_+1);
	y_.resize(N_+1);
	vx_.resize(N_+1);
	vy_.resize(N_+1);
}

void CSolarSystem::add(CPlanet NewPlanet){

	planets_ += 1;
	planet_list_.push_back(NewPlanet);

	// add one column for each planet in the solar system
	for(int i = 0; i < N_+1; i++){
		x[i].resize(planets_);
		y[i].resize(planets_);
		vx[i].resize(planets_);
		vy[i].resize(planets_);
	}
}

// calculates distance between jth and kth planets at ith time step
double CSolarSystem::distance(int i, int j, int k){

	double x, y;

	x = planet_list_[j].x_[i] - planet_list_[k].x_[i];
	y = planet_list_[j].y_[i] - planet_list_[k].y_[i];

	return sqrt(x*x+y*y);
}

// calculates the acceleration in x-direction of the jth planet 
// due to the other planets in the system at the ith time step
double CSolarSystem::ax(int i, int j){

	double m, x, r, ax = 0.0;

	for(int k = 0; k < planets_; j++){
		m = planet_list_[k].m_;
		x = planet_list_[k].x_[i] - planet_list_[j].x_[i];
		r = planet_list_[i].distance(i, j, k);
		ax += 4.0*pi*pi*m*x/(r*r*r);
	}

	return ax;
}

// calculates the acceleration in y-direction of the jth planet 
// due to the other planets in the system at the ith time step
double CSolarSystem::ay(int i, int j){

	double m, x, r, ay = 0.0;

	for(int k = 0; k < planets_; j++){
		m = planet_list_[k].m_;
		y = planet_list_[k].y_[i] - planet_list_[j].y_[i];
		r = planet_list_[i].distance(i, j, k);
		ax += 4.0*pi*pi*m*y/(r*r*r);
	}

	return ax;
}


// sun is fixed at (0,0)
// first planet in list is sun
// solves for trajectory of second planet in list (earth)
void CSolarSystem::solve_euler(string filename){

	vector<double> t(N_+1);
	vector<double> x(N_+1), y(N_+1), vx(N_+1), vy(N_+1);
	double m, v, r, KE, PE;

	// set-up time vector
	t[0] = t0_;
	t[N_] = tf_;
	for(int i = 1; i < N_; i++){ t[i] = t0_+i*h_; }

	// initial conditions of earth
	x[0] = planet_list_[1].x0_[0];
	y[0] = planet_list_[1].x0_[1];
	vx[0] = planet_list_[1].v0_[0];
	vy[0] = planet_list_[1].v0_[1];
	m = planet_list_[1].m_;


	// euler's method
	for(int i = 1; i <= N_; i++){
		x[i] = x[i-1] + h_*vx[i-1];
		y[i] = y[i-1] + h_*vy[i-1];
		vx[i] = vx[i-1] + h_*fx(x[i-1],y[i-1]);
		vy[i] = vy[i-1] + h_*fy(x[i-1],y[i-1]);
	}

	// write to file
	cout << "Writing trajectories to '" << filename << "'... ";
	ofstream ofile;
	ofile.open(filename);
	ofile << "# t0 = " << t0_ << endl;
	ofile << "# tf = " << tf_ << endl;
	ofile << "# N = " << N_ << endl;
	ofile << "# time (year), x (AU), y (AU), kinetic energy, potential energy, total energy" << endl;
	for(int i = 0; i <= N_; i++){
		KE = 0.5*m*(vx[i]*vx[i]+vy[i]*vy[i]);
		PE = -4.0*pi*pi*m/(x[i]*x[i]+y[i]*y[i]);
		ofile << left << setw(14) << setprecision(7) << t[i];
		ofile << left << setw(14) << setprecision(7) << x[i];
		ofile << left << setw(14) << setprecision(7) << y[i];
		ofile << left << setw(14) << setprecision(7) << KE;
		ofile << left << setw(14) << setprecision(7) << PE;
		ofile << left << setw(14) << setprecision(7) << KE+PE << endl;
	}

	// close file
	ofile.close();
	cout << "done!" << endl;

}

// solves for trajectories and energies of all planets
// sun is FIXED at (0,0)
void CSolarSystem::solve_euler(string filename){

	double m, v, r, KE, PE;


	// initial conditions of planets (exclude sun = planet_list_[0])
	for(int j = 1; j < planets_; j++){
		x_[0][j] = planet_list_[j].x0_[0];
		y_[0][j] = planet_list_[j].x0_[1];
		vx_[0][j] = planet_list_[j].v0_[0];
		vy_[0][j] = planet_list_[j].v0_[1];
	}

	// euler's method
	for(int i = 1; i <= N_; i++){
		for(int j = 1; j < planets_; j++){
			x_[i][j] = x_[i-1][j] + h_*vx_[i-1][j];
			y_[i][j] = y_[i-1][j] + h_*vy_[i-1][j];
			vx_[i][j] = vx_[i-1][j] + h_*ax(i-1,j);
			vy_[i][j] = vy_[i-1][j] + h_*ay(i-1,j);
		}
	}

	// write to file
	cout << "Writing trajectories to '" << filename << "'... ";
	ofstream ofile;
	ofile.open(filename);
	ofile << "# t0 = " << t0_ << endl;
	ofile << "# tf = " << tf_ << endl;
	ofile << "# N = " << N_ << endl;
	ofile << "# time (year), x (AU), y (AU), kinetic energy, potential energy, total energy" << endl;
	for(int i = 0; i <= N_; i++){
		KE = 0.5*m*(vx[i]*vx[i]+vy[i]*vy[i]);
		PE = -4.0*pi*pi*m/(x[i]*x[i]+y[i]*y[i]);
		ofile << left << setw(14) << setprecision(7) << t[i];
		ofile << left << setw(14) << setprecision(7) << x[i];
		ofile << left << setw(14) << setprecision(7) << y[i];
		ofile << left << setw(14) << setprecision(7) << KE;
		ofile << left << setw(14) << setprecision(7) << PE;
		ofile << left << setw(14) << setprecision(7) << KE+PE << endl;
	}

	// close file
	ofile.close();
	cout << "done!" << endl;

}

void CSolarSystem::compare_euler(string filename, int maxpower){

	// solves system using euler's method for different number of mesh points
	// writes to files with corresponding power of 10
	for(int i = 2; i <= maxpower; i++){
		N_ = pow(10,i);
		h_ = (tf_-t0_)/N_;
		string outfile = filename + to_string(i) + ".dat";
		this->solve_euler(outfile);
	}
}

// sun is fixed at (0,0)
// first planet in list is sun
// solves for trajectory of second planet in list (earth)
void CSolarSystem::solve_vv(string filename){

	vector<double> t(N_+1);
	vector<double> x(N_+1), y(N_+1), vx(N_+1), vy(N_+1);
	double m, v, r, KE, PE;

	// set-up time vector
	t[0] = t0_;
	t[N_] = tf_;
	for(int i = 1; i < N_; i++){ t[i] = t0_+i*h_; }

	// initial conditions of earth
	x[0] = planet_list_[1].x_[0];
	y[0] = planet_list_[1].x_[1];
	vx[0] = planet_list_[1].v_[0];
	vy[0] = planet_list_[1].v_[1];
	m = planet_list_[1].m_;


	// velocity verlet
	for(int i = 1; i <= N_; i++){
		x[i] = x[i-1] + h_*vx[i-1] + 0.5*h_*h_*fx(x[i-1],y[i-1]);
		y[i] = y[i-1] + h_*vy[i-1] + 0.5*h_*h_*fy(x[i-1],y[i-1]);
		vx[i] = vx[i-1] + 0.5*h_*( fx(x[i],y[i]) + fx(x[i-1],y[i-1]) );
		vy[i] = vy[i-1] + 0.5*h_*( fy(x[i],y[i]) + fy(x[i-1],y[i-1]) );
	}

	// write to file
	cout << "Writing trajectories to '" << filename << "'... ";
	ofstream ofile;
	ofile.open(filename);
	ofile << "# t0 = " << t0_ << endl;
	ofile << "# tf = " << tf_ << endl;
	ofile << "# N = " << N_ << endl;
	ofile << "# time (year), x (AU), y (AU), kinetic energy, potential energy, total energy" << endl;
	for(int i = 0; i <= N_; i++){
		KE = 0.5*m*(vx[i]*vx[i]+vy[i]*vy[i]);
		PE = -4.0*pi*pi*m/(x[i]*x[i]+y[i]*y[i]);
		ofile << left << setw(14) << setprecision(7) << t[i];
		ofile << left << setw(14) << setprecision(7) << x[i];
		ofile << left << setw(14) << setprecision(7) << y[i];
		ofile << left << setw(14) << setprecision(7) << KE;
		ofile << left << setw(14) << setprecision(7) << PE;
		ofile << left << setw(14) << setprecision(7) << KE+PE << endl;
	}

	// close file
	ofile.close();
	cout << "done!" << endl;
}

void CSolarSystem::compare_vv(string filename, int maxpower){

	// solves system using the velocity verlet method for different number of mesh points
	// writes to files with corresponding power of 10
	for(int i = 2; i <= maxpower; i++){
		N_ = pow(10,i);
		h_ = (tf_-t0_)/N_;
		string outfile = filename + to_string(i) + ".dat";
		this->solve_vv(outfile);
	}
}

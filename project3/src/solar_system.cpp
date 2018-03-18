#include "solar_system.h"

CSolarSystem::CSolarSystem(){
	N_ = 1000;
	planets_ = 0;
	t0_ = 0;
	tf_ = 10.0;
	h_ = (tf_-t0_)/N_;
}

CSolarSystem::CSolarSystem(int N, double t0, double tf){
	N_ = N;
	planets_ = 0;
	t0_ = t0;
	tf_ = tf;
	h_ = (tf_-t0_)/N_;	
}

void CSolarSystem::add(CPlanet NewPlanet){
	planets_ += 1;
	planet_list_.push_back(NewPlanet);
}

// sun is fixed at (0,0)
// first planet in list is sun
// solves for trajectory of second planet in list
void CSolarSystem::solve_euler(string filename){

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
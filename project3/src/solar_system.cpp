#include "solar_system.h"

CSolarSystem::CSolarSystem(){

	N_ = 1000;
	dim_ = 2;
	planets_ = 0;
	h_ = 10.0/N_;
	CM_frame_ = false;
	comp_time_ = 0.0;

	// set-up time vector
	t_.resize(N_+1);
	t_[0] = 0.0;
	t_[N_] = 10.0;
	for(int i = 1; i < N_; ++i){ t_[i] = i*h_; }

	// set-up matrices to store data
	x_.resize(N_+1);
	y_.resize(N_+1);
	vx_.resize(N_+1);
	vy_.resize(N_+1);
}

CSolarSystem::CSolarSystem(int N, int dim, double t0, double tf, bool CM_frame){

	N_ = N;
	dim_ = dim;
	planets_ = 0;
	h_ = (tf-t0)/N_;
	CM_frame_ = CM_frame;
	comp_time_ = 0.0;

	// set-up time vector
	t_.resize(N_+1);
	t_[0] = t0;
	t_[N_] = tf;
	for(int i = 1; i < N_; ++i){ t_[i] = t0+i*h_; }

	// set-up matrices to store data
	x_.resize(N_+1);
	y_.resize(N_+1);
	vx_.resize(N_+1);
	vy_.resize(N_+1);

	// extra matrices if 3D
	if(dim_ == 3){
		z_.resize(N_+1);
		vz_.resize(N_+1);
	}
}

void CSolarSystem::add(CPlanet NewPlanet){

	planets_ += 1;
	planet_list_.push_back(NewPlanet);

	// add one column for each body in the solar system
	for(int i = 0; i < N_+1; ++i){

		x_[i].resize(planets_);
		y_[i].resize(planets_);
		vx_[i].resize(planets_);
		vy_[i].resize(planets_);

		if(dim_ == 3){
			z_[i].resize(planets_);
			vz_[i].resize(planets_);
		}
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
	for(int i = 1; i < N_; ++i){ t_[i] = t0+i*h_; }

	// update matrices to store data
	x_.resize(N_+1);
	y_.resize(N_+1);
	vx_.resize(N_+1);
	vy_.resize(N_+1);

	if(dim_ == 3){
		z_.resize(N_+1);
		vz_.resize(N_+1);
	}

	// update number of columns
	for(int i = 0; i <= N_; ++i){

		x_[i].resize(planets_);
		y_[i].resize(planets_);
		vx_[i].resize(planets_);
		vy_[i].resize(planets_);

		if(dim_ == 3){
			z_[i].resize(planets_);
			vz_[i].resize(planets_);
		}
	}
}

// calculates distance between jth and kth planets at ith time step
double CSolarSystem::distance(int i, int j, int k){

	double x, y, z = 0;

	x = x_[i][j]-x_[i][k];
	y = y_[i][j]-y_[i][k];
	if(dim_ == 3) z_[i][j]-z_[i][k];

	return sqrt(x*x+y*y+z*z);
}

// calculates the acceleration of the jth planet due to 
// the other planets in the system at the ith time step
void CSolarSystem::get_acceleration(int i, int j, double& ax, double& ay, double& az){

	double m, r, r3;

	ax = 0.0;
	ay = 0.0;
	az = 0.0;
	for(int k = 0; k < planets_; k++){
		if(k != j) {
			m = planet_list_[k].m_;
			r = distance(i, j, k);
			r3 = r*r*r;
			ax += prefactor*m*(x_[i][k]-x_[i][j])/r3;
			ay += prefactor*m*(y_[i][k]-y_[i][j])/r3;	
			if(dim_ == 3) az += prefactor*m*(z_[i][k]-z_[i][j])/r3;	
		}
	}	
}

// calculates the energy of the jth planet at the ith time step
void CSolarSystem::get_energy(int i, int j, double& KE, double& PE){

	double M, m, r, v2;

	M = planet_list_[j].m_;
	v2 = vx_[i][j]*vx_[i][j]+vy_[i][j]*vy_[i][j];
	if(dim_ == 3) v2 += vz_[i][j]*vz_[i][j]; 
	KE = 0.5*M*v2;

	PE = 0.0;
	for(int k = 0; k < planets_; k++){
		if(k != j){
			m = planet_list_[k].m_;
			r = distance(i, j, k);
			PE += -prefactor*M*m/r;
		}
	}
}

// initial conditions of planets are copied into matrices
// if CM frame is chosen, initial positions are shifted and
// sun is given initial velocity so that total momentum = 0
void CSolarSystem::initialize(){

	// initial conditions of planets (including sun)
	for(int j = 0; j < planets_; j++){

		x_[0][j] = planet_list_[j].x0_[0];
		y_[0][j] = planet_list_[j].x0_[1];
		vx_[0][j] = planet_list_[j].v0_[0];
		vy_[0][j] = planet_list_[j].v0_[1];

		if(dim_ == 3){
			z_[0][j] = planet_list_[j].x0_[2];
			vz_[0][j] = planet_list_[j].v0_[2];		
		}
	}

	// shift initial positions if center-of-mass option is chosen
	// change initial velocity of sun so that total momentum = 0
	CM_[0] = 0.0;
	CM_[1] = 0.0;
	CM_[2] = 0.0;
	if(CM_frame_){

		// find center-of-mass position
		double m, xsum = 0.0, ysum = 0.0, zsum = 0.0;
		for(int j = 0; j < planets_; j++){

			m = planet_list_[j].m_;
			xsum += x_[0][j];
			ysum += y_[0][j];
			CM_[0] += m*x_[0][j];
			CM_[1] += m*y_[0][j];

			if(dim_ == 3){
				zsum += z_[0][j];
				CM_[2] += m*z_[0][j];
			}
		}
		if(xsum != 0) CM_[0] = CM_[0]/xsum;
		if(ysum != 0) CM_[1] = CM_[1]/ysum;
		if(zsum != 0) CM_[2] = CM_[2]/zsum;

		cout << "CM = (" << CM_[0] << ", " << CM_[1] << ", " << CM_[2] << ")\n" << endl;

		// shift
		for(int j = 0; j < planets_; j++){
			x_[0][j] -= CM_[0];
			y_[0][j] -= CM_[1];
			if(dim_ == 3) z_[0][j] -= CM_[2];
		}

		// zero out momentum
		double px = 0.0, py = 0.0, pz = 0.0;
		for(int j = 1; j < planets_; j++){
			m = planet_list_[j].m_;
			px += m*vx_[0][j];
			py += m*vy_[0][j];
			if(dim_ == 3) pz += m*vz_[0][j];
		}
		m = planet_list_[0].m_;
		vx_[0][0] = -px/m;
		vy_[0][0] = -py/m;
		if(dim_ == 3) vz_[0][0] = -pz/m;
	}
}


void CSolarSystem::solve_euler(){

	double a[3];
	int j0;

	// start timer
	clock_t initial, final;
	initial = clock();

	// forward euler
	for(int i = 1; i <= N_; ++i){

		if(CM_frame_) j0 = 0;
		else{
			j0 = 1;

			// fix the sun at (0,0)
			x_[i][0] = 0.0;
			y_[i][0] = 0.0;
			vx_[i][0] = 0.0;
			vy_[i][0] = 0.0;
			if(dim_ == 3){
				z_[i][0] = 0.0;
				vz_[i][0] = 0.0;
			}			
		}

		// calculate orbits
		for(int j = j0; j < planets_; j++){

			x_[i][j] = x_[i-1][j] + h_*vx_[i-1][j];
			y_[i][j] = y_[i-1][j] + h_*vy_[i-1][j];
			if(dim_ == 3) z_[i][j] = z_[i-1][j] + h_*vz_[i-1][j];
			get_acceleration(i-1, j, a[0], a[1], a[2]);
			vx_[i][j] = vx_[i-1][j] + h_*a[0];
			vy_[i][j] = vy_[i-1][j] + h_*a[1];
			if(dim_ == 3) vz_[i][j] = vz_[i-1][j] + h_*a[2];
		}
	}

	final = clock();
	comp_time_ = (final-initial)/((double) CLOCKS_PER_SEC);
}

void CSolarSystem::solve_vv(){

	double a1[3], a2[3];
	int j0;

	// start timer
	clock_t initial, final;
	initial = clock();

	// velocity verlet
	for(int i = 1; i <= N_; ++i){

		if(CM_frame_) j0 = 0;
		else{
			j0 = 1;

			// fix the sun at (0,0)
			x_[i][0] = 0.0;
			y_[i][0] = 0.0;
			vx_[i][0] = 0.0;
			vy_[i][0] = 0.0;	
			if(dim_ == 3){
				z_[i][0] = 0.0;
				vz_[i][0] = 0.0;
			}			
		}

		// calculate orbits
		for(int j = j0; j < planets_; j++){

			get_acceleration(i-1, j, a1[0], a1[1], a1[2]);
			x_[i][j] = x_[i-1][j] + h_*vx_[i-1][j] + 0.5*h_*h_*a1[0];
			y_[i][j] = y_[i-1][j] + h_*vy_[i-1][j] + 0.5*h_*h_*a1[1];
			if(dim_ == 3) z_[i][j] = z_[i-1][j] + h_*vz_[i-1][j] + 0.5*h_*h_*a1[2];

			get_acceleration(i, j, a2[0], a2[1], a2[2]);
			vx_[i][j] = vx_[i-1][j] + 0.5*h_*(a1[0]+a2[0]);
			vy_[i][j] = vy_[i-1][j] + 0.5*h_*(a1[1]+a2[1]);
			if(dim_ == 3) vz_[i][j] = vz_[i-1][j] + 0.5*h_*(a1[2]+a2[2]);
		}	
	}

	final = clock();
	comp_time_ = (final-initial)/((double) CLOCKS_PER_SEC);
}

void CSolarSystem::write_orbits(string filename){

	string Xfile = filename + "_X.dat";
	string Yfile = filename + "_Y.dat";
	string Zfile = filename + "_Z.dat";
	
	// make heading to label columns of file
	string heading = "# time";
	for(int j = 0; j < planets_; j++){
		heading += ", " + planet_list_[j].name_;
	}

	cout << "Writing x-coordinates to '" << Xfile << "'... " << endl;
	cout << "Writing y-coordinates to '" << Yfile << "'... " << endl;
	if(dim_ == 3) cout << "Writing z-coordinates to '" << Zfile << "'... " << endl;

	// print X and Y first
	ofstream outX, outY;
	outX.open(Xfile);
	outY.open(Yfile);

	outX << "# N = " << N_ << endl;
	outX << "# h = " << h_ << endl;
	outX << "# t0 = " << t_[0] << endl;
	outX << "# tf = " << t_[N_] << endl;
	outX << "# computation time = " << comp_time_ << " s" << endl;
	outX << heading << endl;

	outY << "# N = " << N_ << endl;
	outY << "# h = " << h_ << endl;
	outY << "# t0 = " << t_[0] << endl;
	outY << "# tf = " << t_[N_] << endl;
	outY << "# computation time = " << comp_time_ << " s" << endl;
	outY << heading << endl;

	// limit number of points to print
	int istep = 1;
	if(N_ > 1E5) istep = round(N_/1E5);

	// write positions of all planets to file
	for(int i = 0; i <= N_; i += istep){

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

	// now print Z
	if(dim_ == 3){

		ofstream outZ;
		outZ.open(Zfile);

		outZ << "# N = " << N_ << endl;
		outZ << "# h = " << h_ << endl;
		outZ << "# t0 = " << t_[0] << endl;
		outZ << "# tf = " << t_[N_] << endl;
		outZ << "# computation time = " << comp_time_ << " s" << endl;
		outZ << heading << endl;

		// write positions of all planets to file
		for(int i = 0; i <= N_; i += istep){

			outZ << left << setw(14) << setprecision(7) << t_[i];

			// positions of planets
			for(int j = 0; j < planets_; j++){
				outZ << left << setw(14) << setprecision(7) << z_[i][j];
			}
			outZ << endl;
		}
	}
}

void CSolarSystem::write_energies(string filename){

	double KE, PE;
	string Efile = filename + "_E.dat";

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
	outE << "# computation time = " << comp_time_ << " s" << endl;


	outE << heading << endl;

	// limit number of points to print
	int istep = 1;
	if(N_ > 1E5) istep = round(N_/1E5);

	// write all energies to file
	for(int i = 0; i <= N_; i += istep){

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
	initialize();
	for(int n = 3; n <= maxpower; n++){

		cout << "N = 2E" << n << ":" << endl;
		string filename = systemname + "_euler" + to_string(n);
		change_N(2*pow(10,n));
		solve_euler();
		write_orbits(filename);
		write_energies(filename);
		cout << endl;		
	}
}

// writes orbits and energies to files for different N using velocity verlet method
void CSolarSystem::compare_vv(string systemname, int maxpower){

	cout << "\n*** VELOCITY VERLET ***\n" << endl;
	initialize();
	for(int n = 1; n <= maxpower; n++){

		cout << "N = 2E" << n << ":" << endl;
		string filename = systemname + "_vv" + to_string(n);
		change_N(2*pow(10,n));
		solve_vv();
		write_orbits(filename);
		write_energies(filename);
		cout << endl;		
	}
}



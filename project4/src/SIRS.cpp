#include "SIRS.h"

CInfectedPopulation::CInfectedPopulation(){
	N_ = 1E5;
	a_ = 1.0;
	b_ = 1.0;
	c_ = 1.0;
	dt_ = 1E-5;
}

CInfectedPopulation::CInfectedPopulation(int N, double a, double b, double c){
	N_ = N;
	a_ = a;
	b_ = b;
	c_ = c;

	// calculate sufficiently small time step
	dt_ = 1.0;
	if(4.0/(a*N) < dt_) dt_ = 4.0/(a*N);
	if(1.0/(b*N) < dt_) dt_ = 1.0/(b*N);
	if(1.0/(c*N) < dt_) dt_ = 1.0/(c*N);
}

void CInfectedPopulation::deterministic_SIRS(string filename, double S0, double I0, double tf){

	ofstream outfile;
	outfile.open(filename + ".dat");

	cout << "write to ---> " << "'" << filename+".dat'" << endl;

	double S = S0, I = I0, R = N_-S0-I0;
	double S_k1, S_k2, I_k1, I_k2, R_k1, R_k2;

		outfile << "# N = " << N_ << endl;
		outfile << "# (S0, I0, R0) = (" << S << ", " << I << ", " << R << ")" << endl;
		outfile << "# (a, b, c) = (" << a_ << ", " << b_ << ", " << c_ << ")" << endl;

	for(double t = 0.0; t < tf; t += dt_){

		// print 
		outfile << t << "\t" << S << "\t" << I << "\t" << R << endl;

		// use RK2 method to calculate next S, I, R
		S_k1 = dt_*(b_*R-a_*S*I/N_);
		I_k1 = dt_*(a_*S*I/N_-c_*I);
		R_k1 = dt_*(c_*I-b_*R);

		S_k2 = dt_*(b_*(R+0.5*R_k1)-a_*(S+0.5*S_k1)*(I+0.5*I_k1)/N_);
		I_k2 = dt_*(a_*(S+0.5*S_k1)*(I+0.5*I_k1)/N_-c_*(I+0.5*I_k1));
		R_k2 = dt_*(c_*(I+0.5*I_k1)-b_*(R+0.5*R_k1));

		S += S_k2;
		I += I_k2;
		R += R_k2;
	}

	outfile.close();
}

void CInfectedPopulation::generate_phaseportrait(string filename, double tf){

	int points = 4, count = 0;
	double S0, I0, offset = 5.0;
	double spacing = (N_-3.0*offset)/points;

	// vertical wall of triangle
	S0 = offset;
	for(I0 = offset; I0 <= N_-2.0*offset; I0 += spacing){
		deterministic_SIRS(filename+to_string(count), S0, I0, tf);
		count += 1;
	}

	// horizontal wall
	I0 = offset;
	for(S0 = offset+spacing; S0 <= N_-2.0*offset; S0 += spacing){
		deterministic_SIRS(filename+to_string(count), S0, I0, tf);
		count += 1;		
	}

	// hypotenuse
	for(S0 = offset+spacing; S0 <= offset+(points-1)*spacing; S0 += spacing){
		I0 = 4.0*spacing-S0;
		deterministic_SIRS(filename+to_string(count), S0, I0, tf);
		count += 1;		
	}
}

void CInfectedPopulation::montecarlo_SIRS(string filename, int nsamples, int S0, int I0, double tf){

	mt19937 generator;
	uniform_real_distribution<double> distribution(0.0, 1.0);

	vector<double> time;
	for(double t = 0.0; t < tf; t += dt_) time.push_back(t);
	int ntimes = time.size();

	vector<double> avgS(ntimes,0.0), avgI(ntimes,0.0), avgR(ntimes,0.0);

	// create random samples
	for(int n = 0; n < nsamples; ++n){

		ofstream outfile;
		outfile.open(filename + to_string(n) + ".dat");
		cout << "write to ---> " << "'" << filename+to_string(n)+".dat'" << endl;

		int S = S0, I = I0, R = N_-S0-I0;

		outfile << "# N = " << N_ << endl;
		outfile << "# (S0, I0, R0) = (" << S << ", " << I << ", " << R << ")" << endl;
		outfile << "# (a, b, c) = (" << a_ << ", " << b_ << ", " << c_ << ")" << endl;

		for(int i = 0; i < ntimes; ++i){

			// print 
			outfile << time[i] << "\t" << S << "\t" << I << "\t" << R << endl;
			

			// calculate averages
			avgS[i] += S/(double) nsamples;
			avgI[i] += I/(double) nsamples;
			avgR[i] += R/(double) nsamples;

			// keep-or-reject
			if(distribution(generator) < a_*S*I*dt_/N_){ I += 1; S -= 1; }
			if(distribution(generator) < b_*R*dt_){ S += 1; R -= 1; }
			if(distribution(generator) < c_*I*dt_){ R += 1; I -= 1; }
		}
		outfile.close();
	}

	double t, S, I, R;
	string dummy;
	vector<double> sigma(ntimes,0.0);

	// calculate standard deviation
	for(int n = 0; n < nsamples; ++n){

		ifstream infile;
		infile.open(filename + to_string(n) + ".dat");

		// ignore first three lines
		getline(infile, dummy);
		getline(infile, dummy);
		getline(infile, dummy);

		while(!infile.eof()){
			infile >> t >> S >> I >> R >> endl;
			infile >> S;
			infile >> I;
			infile >> R;
			cout << t << "\t" << S << endl;
		}

	

		infile.close();

	}


}



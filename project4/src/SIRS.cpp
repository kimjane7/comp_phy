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



void CInfectedPopulation::deterministic_SIRS(string filename, double S0, double I0, double R0, double tf){

	ofstream outfile;
	outfile.open(filename + ".dat");

	double S = S0, I = I0, R = R0;
	double S_k1, S_k2, I_k1, I_k2, R_k1, R_k2;

	for(double t = 0.0; t < tf; t += dt_){

		// print 
		outfile << left << setw(14) << setprecision(7) << t;
		outfile << left << setw(14) << setprecision(7) << S;
		outfile << left << setw(14) << setprecision(7) << I;
		outfile << left << setw(14) << setprecision(7) << R << endl;

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

void CInfectedPopulation::stochastic_SIRS(string filename, int ntrials, int S0, int I0, int R0, double tf){

	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0, 1.0);
	
	for(int i = 0; i < ntrials; ++i){

		ofstream outfile;
		outfile.open(filename + to_string(i) + ".dat");
		cout << "write to ---> " << "'" << filename << to_string(i) << ".dat'" << endl;

		int S = S0, I = I0, R = R0;

		for(double t = 0.0; t < tf; t += dt_){

			// print 
			outfile << left << setw(14) << setprecision(7) << t;
			outfile << left << setw(14) << setprecision(7) << S;
			outfile << left << setw(14) << setprecision(7) << I;
			outfile << left << setw(14) << setprecision(7) << R << endl;

			// keep-or-reject
			if(distribution(generator) < a_*S*I*dt_/N_){ I += 1; S -= 1; }
			if(distribution(generator) < b_*R*dt_){ S += 1; R -= 1; }
			if(distribution(generator) < c_*I*dt_){ R += 1; I -= 1; }
		}

		outfile.close();

	}

}
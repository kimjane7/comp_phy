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
		outfile << "# time, S, I, R" << endl;

	for(double t = 0.0; t < tf; t += dt_){

		// print 
		outfile << t << "\t" << S << "\t" << I << "\t" << R << endl;

		// use RK2 method to calculate next S, I, R
		S_k1 = dt_*(c_*R-a_*S*I/N_);
		I_k1 = dt_*(a_*S*I/N_-b_*I);
		R_k1 = dt_*(b_*I-c_*R);

		S_k2 = dt_*(c_*(R+0.5*R_k1)-a_*(S+0.5*S_k1)*(I+0.5*I_k1)/N_);
		I_k2 = dt_*(a_*(S+0.5*S_k1)*(I+0.5*I_k1)/N_-b_*(I+0.5*I_k1));
		R_k2 = dt_*(b_*(I+0.5*I_k1)-c_*(R+0.5*R_k1));

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
	uniform_real_distribution<double> rand01(0.0, 1.0);

	vector<double> time;
	for(double t = 0.0; t < tf; t += dt_) time.push_back(t);
	int ntimes = time.size();

	vector<double> avgS(ntimes,0.0), avgI(ntimes,0.0), avgR(ntimes,0.0);

	// create random samples
	for(int n = 0; n < nsamples; ++n){

		ofstream outfile;
		outfile.open(filename+to_string(n)+".dat");
		cout << "write to ---> " << "'" << filename+to_string(n)+".dat'" << endl;

		int S = S0, I = I0, R = N_-S0-I0;

		outfile << "# N = " << N_ << endl;
		outfile << "# (S0, I0, R0) = (" << S << ", " << I << ", " << R << ")" << endl;
		outfile << "# (a, b, c) = (" << a_ << ", " << b_ << ", " << c_ << ")" << endl;
		outfile << "# time, S, I, R" << endl;

		for(int i = 0; i < ntimes; ++i){

			// write to file
			outfile << time[i] << "\t" << S << "\t" << I << "\t" << R << endl;
			
			// calculate averages
			avgS[i] += S/(double) nsamples;
			avgI[i] += I/(double) nsamples;
			avgR[i] += R/(double) nsamples;

			// keep-or-reject
			if(rand01(generator) < a_*S*I*dt_/N_){ I += 1; S -= 1; }
			if(rand01(generator) < b_*I*dt_){ R += 1; I -= 1; }
			if(rand01(generator) < c_*R*dt_){ S += 1; R -= 1; }
		}
		outfile.close();
	}

	double t, S, I, R;
	string dummy;
	vector<double> varS(ntimes,0.0);
	vector<double> varI(ntimes,0.0);
	vector<double> varR(ntimes,0.0);

	// calculate variances
	for(int n = 0; n < nsamples; ++n){

		ifstream infile;
		infile.open(filename + to_string(n) + ".dat");

		// ignore first four lines
		getline(infile, dummy);
		getline(infile, dummy);
		getline(infile, dummy);
		getline(infile, dummy);

		for(int i = 0; i < ntimes; ++i){

			infile >> t >> S >> I >> R;

			varS[i] += (avgS[i]-S)*(avgS[i]-S)/(nsamples-1);
			varI[i] += (avgI[i]-I)*(avgI[i]-I)/(nsamples-1);
			varR[i] += (avgR[i]-R)*(avgR[i]-R)/(nsamples-1);
		}

		infile.close();
	}

	// write averages and standard deviations to file
	ofstream outfile;
	outfile.open(filename+to_string(nsamples)+"_stats.dat");
	cout << "write to ---> " << "'" << filename+to_string(nsamples)+"_stats.dat" << endl;

	outfile << "# N = " << N_ << endl;
	outfile << "# (S0, I0, R0) = (" << S << ", " << I << ", " << R << ")" << endl;
	outfile << "# (a, b, c) = (" << a_ << ", " << b_ << ", " << c_ << ")" << endl;
	outfile << "# nsamples = " << nsamples << endl;
	outfile << "# time, avg S, sigma S, avg I, sigma I, avg R, sigma R" << endl;

	for(int i = 0; i < ntimes; ++i){
		outfile << time[i] << "\t" << avgS[i] << "\t" << sqrt(varS[i]);
		outfile << "\t" << avgI[i] << "\t" << sqrt(varI[i]);
		outfile << "\t" << avgI[i] << "\t" << sqrt(varI[i]) << endl;
	}

	outfile.close();
}

void CInfectedPopulation::lattice_SIRS(string filename, int nsamples, int S0, int I0, double tf){

	mt19937 generator;
	uniform_real_distribution<double> rand01(0.0, 1.0);

	int L = (int) sqrt(N_);
	uniform_int_distribution<int> randint(0, L-1);

	double sum = a_+b_+c_;
	double alpha = a_/sum;    // P(S->I)
	double beta = b_/sum;     // P(I->R)
	double gamma = c_/sum;    // P(R->S)

	vector<double> time;
	for(double t = 0.0; t < tf; t += dt_) time.push_back(t);
	int ntimes = time.size();

	vector<double> avgS(ntimes,0.0), avgI(ntimes,0.0), avgR(ntimes,0.0);

	// create random samples
	for(int n = 0; n < nsamples; ++n){

		ofstream outfile;
		outfile.open(filename+to_string(n)+".dat");
		cout << "write to ---> " << "'" << filename+to_string(n)+".dat'" << endl;

		int S = S0, I = I0, R = N_-S0-I0;

		outfile << "# N = " << N_ << endl;
		outfile << "# (S0, I0, R0) = (" << S << ", " << I << ", " << R << ")" << endl;
		outfile << "# (alpha, beta, gamma) = (" << alpha << ", " << beta << ", " << gamma << ")" << endl;
		outfile << "# time, S, I, R" << endl;

		// fill lattice with susceptibles
		vector<vector<int>> state;
		state.resize(L);
		for(int i = 0; i < L; ++i){
			state[i].resize(L);
			for(int j = 0; j < L; ++j) state[i][j] = 0; 
		}

		// randomly place infected sites & store indices
		int count = 0, i, j;
		vector<vector<int>> infected;
		infected.resize(I0);
		while(count < I0){
			i = randint(generator);
			j = randint(generator);
			if(state[i][j] == 0){
				state[i][j] = 1;
				infected[count].resize(2);
				infected[count][0] = i;
				infected[count][1] = j;
				count += 1;
			}
		}

		// randomly place recovered sites & store indices
		count = 0;
		vector<vector<int>> recovered;
		recovered.resize(R);
		while(count < R){
			i = randint(generator);
			j = randint(generator);	
			if(state[i][j] == 0){
				state[i][j] = 2;
				recovered[count].resize(2);
				recovered[count][0] = i;
				recovered[count][1] = j;
				count += 1;
			}
		}

		for(int t = 0; t < ntimes; ++t){

			// write to file
			outfile << time[t] << "\t" << S << "\t" << I << "\t" << R << endl;
			
			// calculate averages
			avgS[t] += S/(double) nsamples;
			avgI[t] += I/(double) nsamples;
			avgR[t] += R/(double) nsamples;

			if(I > 0){

				// pick an infected site
				uniform_int_distribution<int> randI(0, I-1);
				int r = randI(generator);
				i = infected[r][0];
				j = infected[r][1];

				// recover with probability beta
				if(rand01(generator) < beta){
					R += 1;
					I -= 1;
					state[i][j] = 2;
					recovered.resize(R);
					recovered[R-1][0] = i;
					recovered[R-1][1] = j;
					infected.erase(infected.begin()+r);
				}

				// pick a nearest neighbor with periodic boundary conditions
				int iNN, jNN;
				if(rand01(generator) < 0.5){
					if(rand01(generator) < 0.5){ iNN = (i+1)%L; jNN = j; } 
					else { iNN = (i-1)%L; jNN = j; }
				}
				else{
					if(rand01(generator) < 0.5){ iNN = i; jNN = (j+1)%L; } 
					else { iNN = i; jNN = (j-1)%L; }				
				}

				// infect susceptibles with probability alpha
				if(state[iNN][jNN] == 0){
					if(rand01(generator) < alpha){
						I += 1;
						S -= 1;
						state[iNN][jNN] = 1;
						infected.resize(I);
						infected[I-1][0] = iNN;
						infected[I-1][1] = jNN;
					}
				}				
			}

			if(R > 0){

				// pick a recovered site
				uniform_int_distribution<int> randR(0, R-1);
				r = randR(generator);
				i = recovered[r][0];
				j = recovered[r][1];

				// become susceptible with probability gamma
				if(rand01(generator) < gamma){
					recovered.erase(recovered.begin()+r);
				}
			}



		}
		outfile.close();
	}

	double t, S, I, R;
	string dummy;
	vector<double> varS(ntimes,0.0);
	vector<double> varI(ntimes,0.0);
	vector<double> varR(ntimes,0.0);

	// calculate variances
	for(int n = 0; n < nsamples; ++n){

		ifstream infile;
		infile.open(filename + to_string(n) + ".dat");

		// ignore first four lines
		getline(infile, dummy);
		getline(infile, dummy);
		getline(infile, dummy);
		getline(infile, dummy);

		for(int i = 0; i < ntimes; ++i){

			infile >> t >> S >> I >> R;

			varS[i] += (avgS[t]-S)*(avgS[t]-S)/(nsamples-1);
			varI[i] += (avgI[t]-I)*(avgI[t]-I)/(nsamples-1);
			varR[i] += (avgR[t]-R)*(avgR[t]-R)/(nsamples-1);
		}

		infile.close();
	}

	// write averages and standard deviations to file
	ofstream outfile;
	outfile.open(filename+to_string(nsamples)+"_stats.dat");
	cout << "write to ---> " << "'" << filename+to_string(nsamples)+"_stats.dat" << endl;

	outfile << "# N = " << N_ << endl;
	outfile << "# (S0, I0, R0) = (" << S << ", " << I << ", " << R << ")" << endl;
	outfile << "# (a, b, c) = (" << a_ << ", " << b_ << ", " << c_ << ")" << endl;
	outfile << "# nsamples = " << nsamples << endl;
	outfile << "# time, avg S, sigma S, avg I, sigma I, avg R, sigma R" << endl;

	for(int i = 0; i < ntimes; ++i){
		outfile << time[i] << "\t" << avgS[i] << "\t" << sqrt(varS[i]);
		outfile << "\t" << avgI[i] << "\t" << sqrt(varI[i]);
		outfile << "\t" << avgI[i] << "\t" << sqrt(varI[i]) << endl;
	}

	outfile.close();


}

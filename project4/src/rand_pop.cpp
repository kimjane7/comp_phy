#include "rand_pop.h"

CRandPop::CRandPop(){

	N_ = 1000;
	S_ = N_-1;
	I_ = 1;
	R_ = 0;

	L_ = 10.0;
	H_ = 10.0;

}

CRandPop::CRandPop(int N, int I, double L, double H){

	N_ = N;
	S_ = N-I;
	I_ = I;
	R_ = 0;

	L_ = L;
	H_ = H;

	populate();
}

double CRandPop::random01(mt19937 generator){

    static uniform_01<mt19937> dist(generator);
    return dist();
}

void CRandPop::populate(){

	mt19937 generator(time(0));

}
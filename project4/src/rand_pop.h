#ifndef RAND_POP_H
#define RAND_POP_H

#include <string> 
#include <iomanip>
#include <vector>
#include <boost/random.hpp>

using namespace std;
using boost::uniform_01;
using boost::mt19937;

class CRandPop{
public:

	int N_, S_, I_, R_;
	double L_, H_;

	CRandPop();
	CRandPop(int N, int I, double L, double H);
	~CRandPop(){};

	void populate();

};





#endif
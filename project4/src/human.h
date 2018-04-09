#ifndef HUMAN_H
#define HUMAN_H

#include <string> 
#include <iomanip>
#include <vector>

using namespace std;

class CHuman{
public:

	int state;
	double x_, y_;

	CHuman();
	CHuman(int state, double x, double y);
	~CHuman(){};


};





#endif
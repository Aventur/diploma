#include <iostream>
using namespace std;

#ifndef BSVGENERATOR
#define BSVGENERATOR

class LinRand
{
	unsigned __int64 a;
	unsigned __int64 c;
	unsigned __int64 x0;
	unsigned __int64 M;

public:
	LinRand();
	LinRand(unsigned __int64 M, unsigned __int64 a,
		unsigned __int64 c, unsigned __int64 x0);
	~LinRand();
	double next();
	double *next(int N);
	void setx0(__int64 x);
};

class BSV
{
	LinRand *g1;
	LinRand *g2;
	unsigned int K;
	double *table;

public:

	BSV();
	BSV(unsigned int K);
	~BSV();
	void init();
	double next();
	double *next(int N);
	void setx0(__int64 x);
};


#endif
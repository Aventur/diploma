#include "bsv.h"

LinRand::LinRand()
{
	a = c = x0 = M = 0;
}

LinRand::LinRand(unsigned __int64 M, unsigned __int64 a,
		unsigned __int64 c, unsigned __int64 x0)
{
	this->M = M;
	this->a = a;
	this->c = c;
	this->x0 = x0;
}

double LinRand::next()
{
	x0 = (a*x0 + c) % M;
	return (double(x0)/double(M));
}

double* LinRand::next(int N)
{
	double *r = new double[N];
	for (int i = 0; i < N; i++)
		r[i] = this->next();
	return r;
}

LinRand::~LinRand()
{ }

void LinRand::setx0(__int64 x)
{
	this->x0 = x;
}

BSV::BSV()
{
	K = 128;
	init();
}

BSV::BSV(unsigned int K)
{
	this->K = K;
}

void BSV::init()
{
	unsigned __int64 M = INT_MAX + 1;
	g1 = new LinRand(M, 16385, 16777215, 571824);
	g2 = new LinRand(524287, 1048575, 264202, 2);
	table = g1->next(K);
}

double BSV::next()
{
	int s = int(floor(g2->next()*K));
	double r = table[s];
	table[s] = g1->next();
	return r;
}

double* BSV::next(int N)
{
	double *r = new double[N];
	for (int i = 0; i < N; i++)
		r[i] = this->next();
	return r;
}

BSV::~BSV()
{
	if (g1 != NULL) delete g1;
	if (g2 != NULL) delete g2;
}

void BSV::setx0(__int64 x)
{
	g1->setx0(x);
	if (table) delete [] table;
	table = g1->next(K);
}
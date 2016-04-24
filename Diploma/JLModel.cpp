#include "JLModel.h"

JacobsLewisModel::JacobsLewisModel()
{
	ro = 0;
	lambda = NULL;
	P = NULL;
}

JacobsLewisModel::JacobsLewisModel(char *fname)
{
	fstream is;
	is.open(fname, ios::in);
	is >> s;
	is >> L;
	is.ignore(1);
	alphabet = new char[L + 1];
	is.read(alphabet, L);
	alphabet[L] = '\0';
	P = new double[L];
	lambda = new double[s];
	for (int i = 0; i < L; i++)
		is >> P[i];
	for (int i = 0; i < s; i++)
		is >> lambda[i];
	is >> ro;
	buffer = new int[s];
	is.close();
}

void JacobsLewisModel::printModel(ostream &os)
{
	MarkovChainModel::printModel(os);
	printArray(P, L, os);
	printArray(lambda, s, os);
	os << ro;
}

JacobsLewisModel::~JacobsLewisModel()
{
	if (lambda != NULL)
		delete[]lambda;
	if (P != NULL)
		delete[]P;
}

int JacobsLewisModel::nextInitialState()
{
	return nextDiscreteRand(P, L);
}

int JacobsLewisModel::nextState()
{
	int res;
	if (bsv->next()<ro)
		res = buffer[s - nextDiscreteRand(lambda, s) - 1];
	else
		res = nextInitialState();
	shiftBuffer(res);
	return res;
}
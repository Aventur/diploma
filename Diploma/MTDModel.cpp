#include "MTDModel.h"

MTDModel::MTDModel()
{
	Q = NULL;
	lambda = NULL;
	P = NULL;
}

MTDModel::MTDModel(char *fname)
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
	Q = new double*[L];
	for (int i = 0; i < L; i++)
		Q[i] = new double[L];
	for (int i = 0; i < L; i++)
		for (int j = 0; j < L; j++)
			is >> Q[i][j];
	buffer = new int[s];
	curDist = new double[L];
	is.close();
}

MTDModel::~MTDModel()
{
	if (lambda != NULL)
		delete[]lambda;
	if (P != NULL)
		delete[]P;
	if (Q != NULL)
	{
		for (int i = 0; i < L; i++)
			if (Q[i] != NULL) delete[]Q[i];
		delete[]Q;
	}
}

void MTDModel::printModel(ostream &os)
{
	MarkovChainModel::printModel(os);
	printArray(P, L, os);
	printArray(lambda, s, os);
	for (int i = 0; i < L; i++)
		printArray(Q[i], L, os);
}

int MTDModel::nextInitialState()
{
	return nextDiscreteRand(P, L);
}

int MTDModel::nextState()
{
	for (int i = 0; i < L; i++)
	{
		curDist[i] = 0;
		for (int j = 0; j < s; j++)
			curDist[i] += lambda[j] * Q[buffer[j]][i];
	}
	int res = nextDiscreteRand(curDist, L);
	shiftBuffer(res);
	return res;
}
#include "MCModel.h"

MarkovChainModel::MarkovChainModel()
{
	bsv = new BSV();
	bsv->setx0(time(NULL));
	s = 0;
	L = 0;
	alphabet = NULL;
	buffer = NULL;
}

MarkovChainModel::~MarkovChainModel()
{
	if (bsv != NULL)
		delete bsv;
	if (alphabet != NULL)
		delete[]alphabet;
	if (buffer != NULL)
		delete[]buffer;
}

void MarkovChainModel::shiftBuffer(int state)
{
	for (int i = 0; i < (s - 1); i++)
		buffer[i] = buffer[i + 1];
	buffer[s - 1] = state;
}

int MarkovChainModel::indexof(char t)
{
	int i = 0;
	while ((alphabet[i] != t) && (i < L))
		i++;
	if (i == L)
		return (-1);
	return i;
}

void MarkovChainModel::setParam(int p1, int p2)
{}

void MarkovChainModel::estimateModel(istream &is)
{}

void MarkovChainModel::printArray(double *arr, int n, ostream &os)
{
	for (int i = 0; i < n; i++)
		os << arr[i] << " ";
	os << endl;
}

void MarkovChainModel::printArray(int *arr, int n, ostream &os)
{
	for (int i = 0; i < n; i++)
		os << arr[i] << " ";
	os << endl;
}

void MarkovChainModel::printModel(ostream &os, short full)
{
	if (full <= 2)
		os << s << endl;
	os << L << endl;
	os << alphabet << endl;
}

int MarkovChainModel::nextDiscreteRand(double *dist, int n)
{
	int res = 0;
	double a = bsv->next();
	while (res < (n - 1) && a > dist[res])
	{
		a -= dist[res];
		res++;
	}
	return res;
}

void MarkovChainModel::generateSequence(int n, ostream &os)
{
	for (int i = 0; i < s; i++)
	{
		buffer[i] = nextInitialState();
		os << alphabet[buffer[i]];
	}
	for (int i = 0; i < n; i++)
		os << alphabet[nextState()];
}
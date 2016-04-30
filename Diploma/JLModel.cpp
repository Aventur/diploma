#include "JLModel.h"

JacobsLewisModel::JacobsLewisModel()
{
	ro = 0;
	lambda = NULL;
	P = NULL;
	nu = NULL;
}

JacobsLewisModel::JacobsLewisModel(int s, int L, char* alpha)
{
	this->s = s;
	this->L = L;
	this->alphabet = new char[L+1];
	strcpy_s(this->alphabet, L+1, alpha);
	this->alphabet[L] = '\0';
	this->buffer = new int[L];

	ro = .5;
	lambda = new double[s];
	P = new double[L];

	for (int i = 0; i < s; i++)
		lambda[i] = 1. / s;
	for (int i = 0; i < L; i++)
		P[i] = 1. / L;

	nu = NULL;
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
	nu = NULL;
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
	if (nu != NULL)
		delete[]nu;
}

void JacobsLewisModel::estimateNu(istream * is)
{
	if (nu != NULL) delete[] nu;				// заполнение матрицы частот нулями
	nu = new unsigned int[pow(L, s + 1)];
	ZeroMemory(nu, sizeof(int) * pow(L, s + 1));

	Buffer *b = new Buffer(L, s + 1);
	for (int i = 0; i <= s + 1; i++)
		b->shiftBuffer(indexof(is->get()));

	do
	{
		nu[b->getNumber()]++;
		b->shiftBuffer(indexof(is->get()));
	} while (!is->eof());

	delete b;
}

double * JacobsLewisModel::estimateQ(double *Q)
{
	double sum = nu[0];
	int i, j;
	for (i = 1; i < pow(L, s + 1); i++)
	{
		sum += nu[i];
		if (i % L == L - 1)
		{
			for (j = 0; j < L; j++)
				Q[i - j] = nu[i - j] / sum;
			sum = 0;
		}
	}
	return Q;
}

double JacobsLewisModel::estimateInitialRo(double * Q)
{
	Buffer b(L, s + 1);
	double dividend = 0, divisor = 0;
	for (int i = 0; i < L; i++)
	{
		next2(b, (i + 1) % L);
		b[-1] = i;
		for (int j = 0; j < pow(L - 1, s); j++)
			dividend += P[i] * (P[i] - Q[next1(b, i)]);

		dividend -= (1 - P[i]) * (P[i] - Q[next2(b, i)]);
		divisor += pow(L - 1, s) * P[i] * P[i]; + pow(1 - P[i], 2);
	}
	ro = dividend / divisor;
	return ro;
}

double * JacobsLewisModel::estimateInitialLambda(double * Q)
{
	Buffer b(L, s);
	double *C = new double[s];

	double c = 0;					// вычисление C
	for (int i = 0; i < L; i++)
		c += P[i] * P[i];
	c *= 2 * pow(L, s - 1) * ro * ro;

	do
	{
		for (int j = 0; j < s; j++)
			C[j] -= P[b[j]] - Q[L*b.getNumber() + j];
		b++;
	} while (b.getNumber());
	for (int i = 0; i < s; i++)
		C[i] = c + 2 * ro * C[i];

	for (int i = 0; i < s; i++)		// вычисление lambda
	{
		double prod = 0;
		for (int j = 0; j < s; j++)
			prod += (s * (i == j) - 1) * C[j];
		lambda[i] = prod / (2 * s * (L - 1) * pow(L, s - 1) * ro * ro) - 1. / s;
	}

	return lambda;
}

__int64 JacobsLewisModel::next1(Buffer & b, int i)
{
	int j = b.len - 2;
	while (j >= 0)
	{
		if (b[j] < L - 1)
		{
			b[j]++;
			if (b[j] == i) continue;
			break;
		}
		if (!i) b[j] = 1;
		else b[j] = 0;
		j--;
	}
	return b.getNumber();
}

__int64 JacobsLewisModel::next2(Buffer & b, int i)
{
	for (int j = 0; j < b.len; j++)
		b[j] = i;
	return b.getNumber();
}

void JacobsLewisModel::estimateP(istream * is)
{
	unsigned int *frequencies = new unsigned int[L];
	ZeroMemory(frequencies, sizeof(unsigned int) * L);
	
	for (char c = is->get(); !is->eof(); c = is->get())
		frequencies[indexof(c)]++;

	unsigned int sum = frequencies[0];
	for (int i = 1; i < L; i++)
		sum += frequencies[i];

	for (int i = 0; i < L; i++)
		P[i] = double(frequencies[i]) / sum;

	//if (frequencies != NULL) delete[]frequencies;
}

void JacobsLewisModel::estimateInitialParameters(istream * is)
{
	if (nu == NULL) estimateNu(is);
	refresh_string_stream(is);

	// вычислим оценку матрицы Q
	double *Q = new double[pow(L, s + 1)];
	estimateQ(Q);

	estimateP(is);				// оценка начального значения P
	refresh_string_stream(is);

	estimateInitialRo(Q);		// оценка начального значения Ro
	estimateInitialLambda(Q);	// оценка начального значения Lambda

	/*
	next2(b, 1);
	b[-1] = 0;
	for (int i = 0; i < pow(L-1, s) * 4; i++)
	{
		next1(b, 0);
		cout << i << ' ' << b.getNumber() << " ";
		for (int j = 0; j < b.len; j++)
			cout << b[j];
		cout << '\n';
	}*/


	if (Q != NULL) delete[]Q;
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

double JacobsLewisModel::likelihood()
{// TODO!!!
	return 0;
}
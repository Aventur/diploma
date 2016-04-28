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

void JacobsLewisModel::estimateNu(istream & is)
{
	istringstream *sis = get_string_stream(is, alphabet);

	if (nu != NULL) delete[] nu;				// заполнение матрицы частот нулями
	nu = new unsigned int[pow(L, s + 1)];
	ZeroMemory(nu, sizeof(int) * pow(L, s + 1));

	Buffer *b = new Buffer(L, s + 1);
	for (int i = 0; i <= s + 1; i++)
		b->shiftBuffer(indexof(sis->get()));

	do
	{
		nu[b->getNumber()]++;
		b->shiftBuffer(indexof(sis->get()));
	} while (!sis->eof());

	delete sis;
	delete b;
}

void JacobsLewisModel::estimateP(istream & is)
{
	istringstream *sis = get_string_stream(is, alphabet);

	unsigned int *frequencies = new unsigned int[L];
	ZeroMemory(frequencies, sizeof(unsigned int) * L);
	
	for (char c = sis->get(); !sis->eof(); c = sis->get())
		frequencies[indexof(c)]++;

	unsigned int sum = frequencies[0];
	for (int i = 1; i < L; i++)
		sum += frequencies[i];

	for (int i = 0; i < L; i++)
		P[i] = double(frequencies[i]) / sum;

	//if (frequencies != NULL) delete[]frequencies;
	delete sis;
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
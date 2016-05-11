#include "CMModel.h"

void ClassicMarkovModel::estimateNu(istream * is)
{
	if (nu == NULL) nu = new unsigned int[pow(L, s + 1)];
	ZeroMemory(nu, sizeof(unsigned int) * pow(L, s + 1));

	Buffer *b = new Buffer(L, s + 1);
	for (int i = 0; i <= s + 1; i++)
		b->shiftBuffer(indexof(is->get()));

	do
	{
		nu[b->getNumber()]++;
		b->shiftBuffer(indexof(is->get()));
	} while (!is->eof());
	
	refresh_string_stream(is);
	delete b;
}

void ClassicMarkovModel::estimateNu(istream * is, int s_pos, int e_pos)
{
	if (nu == NULL) nu = new unsigned int[pow(L, s + 1)];
	ZeroMemory(nu, sizeof(unsigned int) * pow(L, s + 1));
	is->seekg(s_pos, ios::beg);

	Buffer *b = new Buffer(L, s + 1);
	int pos;
	for (pos = s_pos; pos <= s_pos + s + 1; pos++)
		b->shiftBuffer(indexof(is->get()));

	for (; pos < e_pos; pos++)
	{
		nu[b->getNumber()]++;
		b->shiftBuffer(indexof(is->get()));
	}

	refresh_string_stream(is);
	delete b;
}

double ClassicMarkovModel::get_buffer_likelihood(Buffer & b)
{
	return Q[b.getNumber()];
}

ClassicMarkovModel::ClassicMarkovModel()
{
	s = L = 0;
	alphabet = NULL;
	nu = NULL;
	Q = NULL;
}

ClassicMarkovModel::ClassicMarkovModel(int s, int L, char * alpha)
{
	this->s = s;
	this->L = L;
	this->alphabet = new char[L + 1];
	strcpy_s(this->alphabet, L + 1, alpha);
	this->alphabet[L] = '\0';

	buffer = new int[s];
	nu = new unsigned int[pow(L, s + 1)];
	Q = new double[pow(L, s + 1)];
}

ClassicMarkovModel::ClassicMarkovModel(char * fname)
{
	fstream is;
	is.open(fname, ios::in);
	is >> s;
	is >> L;
	is.ignore(1);
	alphabet = new char[L + 1];
	is.read(alphabet, L);
	alphabet[L] = '\0';
	Q = new double[pow(L, s+1)];
	buffer = new int[s];
	is.close();
	nu = NULL;
}

ClassicMarkovModel::~ClassicMarkovModel()
{
	if (nu != NULL) delete[] nu;
	if (Q != NULL) delete[] Q;
}

double ClassicMarkovModel::estimateQ(istream *is)
{
	estimateNu(is);
	double sum = nu[0];
	int i, j;
	for (i = 1; i < pow(L, s + 1); i++)
	{
		sum += nu[i];
		if (i % L == L - 1)
		{
			for (j = 0; j < L; j++)
				Q[i - j] = sum ? (nu[i - j] / sum) : 1. / L;
			sum = 0;
		}
	}
	return likelihood(is);
}

double ClassicMarkovModel::estimateQ(istream * is, int s_pos, int e_pos)
{
	estimateNu(is, s_pos, e_pos);
	double sum = nu[0];
	int i, j;
	for (i = 1; i < pow(L, s + 1); i++)
	{
		sum += nu[i];
		if (i % L == L - 1)
		{
			for (j = 0; j < L; j++)
				Q[i - j] = sum ? (nu[i - j] / sum) : 1. / L;
			sum = 0;
		}
	}
	return likelihood(is, s_pos, e_pos);
}

double ClassicMarkovModel::likelihood(istream * is)
{
	double lh = 0;
	int cur, i;
	Buffer b(L, s + 1);

	for (i = 0; i < s; i++)
	{
		cur = indexof(is->get());
		b.shiftBuffer(cur);
	}

	for (char c = is->get(); !is->eof(); c = is->get())
	{
		cur = indexof(c);
		if (cur == -1)
			continue;
		b.shiftBuffer(cur);

		lh += log(get_buffer_likelihood(b));
	}

	refresh_string_stream(is);
	return lh;
}

double ClassicMarkovModel::likelihood(istream * is, int s_pos, int e_pos)
{
	double lh = 0;
	int cur, pos, i;
	Buffer b(L, s);
	is->seekg(s_pos, ios::beg);

	for (pos = s_pos; pos < s_pos + s; pos++)
	{
		cur = indexof(is->get());
		b.shiftBuffer(cur);
	}

	for (char c = is->get(); pos < e_pos; c = is->get())
	{
		cur = indexof(c);
		if (cur == -1)
			continue;
		b.shiftBuffer(cur);

		lh += log(get_buffer_likelihood(b));
		pos++;
	}

	refresh_string_stream(is);
	return lh;
}

__int64 ClassicMarkovModel::numberOfParams()
{
	return pow(L, s) * (L - 1);
}

int ClassicMarkovModel::nextInitialState()
{
	double a = bsv->next();
	for (int i = 1; i <= L; i++)
		if (a <= double(i) / L)
			return i - 1;
	return L-1;
}

int ClassicMarkovModel::nextState()
{
	return nextDiscreteRand(&Q[number() * L], L);
}

void ClassicMarkovModel::printModel(ostream & os)
{
	MarkovChainModel::printModel(os);
	for (int i = 0; i < pow(L, s); i++)
	printArray(&Q[i * L], L, os);
	os << endl;
}

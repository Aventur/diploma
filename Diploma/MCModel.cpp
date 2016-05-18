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

int MarkovChainModel::number(char *t)
{
	int r = 0;
	int i = 0;
	while (t[i] != '\0')
	{
		r *= L;
		r += indexof(t[i]);
		i++;
	}
	return r;
}

int MarkovChainModel::number()
{
	int r = 0;
	int i = 0;
	for (i = 0; i < s; r = r * L + buffer[i++]);
	return r;
}

double MarkovChainModel::AIC(istream * is, int s_pos, int e_pos)
{
	return -2. * likelihood(is, s_pos, e_pos) + 2. * (double)numberOfParams();
}

double MarkovChainModel::BIC(istream * is, int s_pos, int e_pos)
{
	int len = e_pos - s_pos;
	return -2. * likelihood(is, s_pos, e_pos) + numberOfParams() * log(len);
}

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

void MarkovChainModel::generateRandomDistribution(double * D, int size)
{
	int i;
	for (i = 0; i < size - 1; i++)
		D[i] = bsv->next();
	D[size - 1] = 1.;

	qsort(D, size, sizeof(double), cmp<double>);
	for (i = size - 1; i > 0; i--)
		D[i] -= D[i - 1];
}

void MarkovChainModel::generateSequence(int n, ostream &os)
{
	for (int i = 0; i < s; i++)
	{
		buffer[i] = nextInitialState();
		os << alphabet[buffer[i]];
	}
	for (int i = 0; i < n - s; i++)
		os << alphabet[nextState()];
}

Buffer::Buffer(int L, int len)
{
	this->L = L;
	this->len = len;
	str = new int[len];
	ZeroMemory(str, sizeof(int) * len);

	beg = 0;
	end = len - 1;
}

Buffer::~Buffer()
{
	if (str != NULL)
		delete[] str;
}

__int64 Buffer::getNumber()
{
	__int64 r = 0;
	for (int i = beg; true; i = (i + 1) % len)
	{
		r = r * L + str[i];
		if (i == end) break;
	}
	return r;
}

void Buffer::shiftBuffer(int i)
{
	str[beg] = i;
	end = beg;
	beg = (beg + 1) % len;
}

int & Buffer::operator[](const int index)
{
	if (index >= len || index < 0) 
		return str[end];
	return str[(beg + index) % len];
}

void Buffer::operator++(int)
{
	int i = len-1;
	while (i >= 0)
	{
		if ((*this)[i] < L - 1)
		{
			(*this)[i]++;
			break;
		}
		(*this)[i] = 0;
		i--;
	}
}



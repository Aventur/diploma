#include "PCModel.h"

PartialConnectionsModel::PartialConnectionsModel()
{
	P = NULL;
	Q = NULL;
	r = 0;
	m = NULL;
}

PartialConnectionsModel::PartialConnectionsModel(char *fname, short full)
{
	fstream is;
	is.open(fname, ios::in);
	if (full <= 2)
		is >> s;
	is >> L;
	is.ignore(1);
	alphabet = new char[L + 1];
	is.read(alphabet, L);
	alphabet[L] = '\0';
	P = new double[L];
	if (full == 0)
		for (int i = 0; i < L; i++)
			is >> P[i];
	if (full <= 2)
	{
		is >> r;
		m = new int[r];
		for (int i = 0; i < r; i++)
			is >> m[i];
		c = pow(double(L), r);
		Q = new double*[c];
		for (int i = 0; i < c; i++)
			Q[i] = new double[L];
	}
	if (full == 0)
		for (int i = 0; i < c; i++)
			for (int j = 0; j < L; j++)
				is >> Q[i][j];
	if (full <= 2)
		buffer = new int[s];
	is.close();
}

PartialConnectionsModel::~PartialConnectionsModel()
{
	if (P != NULL)
		delete[]P;
	if (m != NULL)
		delete[]m;
	if (Q != NULL)
	{
		for (int i = 0; i < c; i++)
			if (Q[i] != NULL) delete[]Q[i];
		delete[]Q;
	}
}

int PartialConnectionsModel::number()
{
	int i = 0;
	for (int j = 0; j < r; j++)
	{
		i *= L;
		i += buffer[m[j]];
	}
	return i;
}

int PartialConnectionsModel::nextInitialState()
{
	return nextDiscreteRand(P, L);
}

int PartialConnectionsModel::nextState()
{
	int i = number();
	int res = nextDiscreteRand(Q[i], L);
	shiftBuffer(res);
	return res;
}

void PartialConnectionsModel::printModel(ostream &os, short full)
{
	MarkovChainModel::printModel(os, full);
	if (full == 0) printArray(P, L, os);
	if (full <= 2)
	{
		os << r << endl;
		if (full <= 1)
		{
			for (int i = 0; i < r; i++)
				os << m[i] << " ";
			os << endl;
		}
	}
	if (full == 0)
		for (int i = 0; i < c; i++)
			printArray(Q[i], L, os);
}

bool PartialConnectionsModel::nextM(int level)
{
	bool res = false;
	if (m[r - level] != s - level)
	{
		m[r - level]++;
		return false;
	}
	else
		if (level != r)
			res = nextM(level + 1);
		else
			return true;
	if (!res)
		m[r - level] = m[r - level - 1] + 1;
	return res;
}

void PartialConnectionsModel::estimateModel(istream &is)
{
	is.seekg(0, ios::end);
	int g = is.tellg();
	is.seekg(0, ios::beg);
	char *file = new char[g];
	is.read(file, g);

	istringstream sis(file);

	int i;
	for (i = 0; i < r; i++)			// начальный вектор M
		m[i] = i;
	int mini = 0;
	double minH = estimateModelQ(sis);
	double H = 0;
	i = 0;
	sis.clear();
	sis.seekg(0, ios::beg);

	while (!nextM(1))				// перебор всех возм. M
	{
		i++;
		H = estimateModelQ(sis);
		if (H < minH)
		{
			minH = H;
			mini = i;
		}
		sis.clear();
		sis.seekg(0, ios::beg);
	}

	for (i = 0; i < r; i++)			// окончательная оценка Q для лучшего M
		m[i] = i;
	for (i = 0; i < mini; i++)
		nextM(1);
	estimateModelQ(sis);

	if (file != NULL) delete[]file;
}

double PartialConnectionsModel::estimateModelQ(istream &is)
{
	if (is.bad()) cout << "Error! Bad input file for estimating PCModel\n";
	for (int i = 0; i < c; i++)		// заполнение Q нулями
		for (int j = 0; j < L; j++)
			Q[i][j] = 0;
	for (int i = 0; i < L; i++)		// заолнение P нулями
		P[i] = 0;

	char t;
	int l = 0;
	int i = 0;
	double H = 0;					// для определения Mr
	double m = 0;
	double sum = 0;

	l = 0;							// начальная подпоследовательность
	while ((l <= s) && (!is.bad()))
	{
		t = is.get();
		i = indexof(t);
		if (i != -1)
		{
			P[i] += 1;
			shiftBuffer(i);
			l++;
		}
	}

	while (true)					// заполнение Q, P частотами
	{
		t = is.get();
		if (is.eof()) break;
		while ((indexof(t) == -1) && (!is.eof()))
			t = is.get();
		if (is.eof()) break;

		i = indexof(t);
		P[i] += 1;
		Q[number()][i] += 1;
		shiftBuffer(i);
	}

	for (i = 0; i < c; i++)			// нормировка Q, P
	{
		sum = 0;
		for (l = 0; l < L; l++)
			sum += Q[i][l];
		if (sum != 0)
			for (l = 0; l < L; l++)
			{
				if (Q[i][l] != 0)
				{
					m = Q[i][l];
					Q[i][l] /= sum;
					H -= m*log(Q[i][l]);
				}
			}
		else
		{
			sum = 1. / L;
			for (l = 0; l < L; l++)
				Q[i][l] = sum;
		}
	}

	sum = 0;
	for (i = 0; i < L; i++)
		sum += P[i];
	for (i = 0; i < L; i++)
		P[i] /= sum;
	H /= sum - s;
	return H;
}

void PartialConnectionsModel::setParam(int p1, int p2)
{
	if (m != NULL)
		delete[]m;
	if (Q != NULL)
	{
		for (int i = 0; i < c; i++)
			if (Q[i] != NULL) delete[]Q[i];
		delete[]Q;
	}
	if (buffer != NULL) delete[]buffer;
	s = p1;
	r = p2;

	m = new int[r];
	c = pow(double(L), r);
	Q = new double*[c];
	for (int i = 0; i < c; i++)
		Q[i] = new double[L];
	buffer = new int[s];
}
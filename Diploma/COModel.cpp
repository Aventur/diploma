#include "COModel.h"

ConditionalOrderModel::ConditionalOrderModel()
{
	P = NULL;
	B = M = K = 0;
	mk = bk = NULL;
	Q = NULL;
}

ConditionalOrderModel::ConditionalOrderModel(char *fname, short full)
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
		is >> B;
	K = pow(double(L), B);
	if (full <= 1)
		is >> M;
	else
		M = K;

	mk = new int[K];
	if (full <= 1)
		for (int i = 0; i < K; i++)
			is >> mk[i];
	else
		for (int i = 0; i < K; i++)
			mk[i] = i;
	bk = new int[K];
	if (full <= 1)
		for (int i = 0; i < K; i++)
			is >> bk[i];
	else
		bk[0] = -1;

	Q = new double*[M*L];
	for (int i = 0; i < M*L; i++)
		Q[i] = new double[L];
	if (full == 0)
		for (int i = 0; i < M*L; i++)
			for (int j = 0; j < L; j++)
				is >> Q[i][j];

	buffer = new int[s];
	is.close();
}

ConditionalOrderModel::~ConditionalOrderModel()
{
	if (P != NULL)
		delete[]P;
	if (mk != NULL)
		delete[]mk;
	if (bk != NULL)
		delete[]bk;
	if (Q != NULL)
	{
		for (int i = 0; i < M*L; i++)
			if (Q[i] != NULL) delete[]Q[i];
		delete[]Q;
	}
}

int ConditionalOrderModel::number()
{
	int i = 0;
	for (int j = s - B; j < s; j++)
	{
		i *= L;
		i += buffer[j];
	}
	return i;
}

void ConditionalOrderModel::setParam(int p1, int p2)
{
	if (mk != NULL)
		delete[]mk;
	if (bk != NULL)
		delete[]bk;
	if (Q != NULL)
	{
		for (int i = 0; i < M*L; i++)
			if (Q[i] != NULL) delete[]Q[i];
		delete[]Q;
	}
	if (buffer != NULL) delete[]buffer;
	s = p1;
	B = p2;

	K = pow(double(L), B);
	M = K;

	mk = new int[K];
	for (int i = 0; i < K; i++)
		mk[i] = i;
	bk = new int[K];
	bk[0] = -1;

	Q = new double*[M*L];
	for (int i = 0; i < M*L; i++)
		Q[i] = new double[L];
	buffer = new int[s];
}

void ConditionalOrderModel::estimateModel(istream &is)
{
	if (bk[0] == -1)
		estimateModelB(is);
	else
		estimateModelQ(is);
}

void ConditionalOrderModel::estimateModelB(istream &is)
{
	if (is.bad()) cout << "Error! Bad input file for estimating PCModel\n";
	for (int i = 0; i < M*L; i++)		// заполнение Q нулями
		for (int j = 0; j < L; j++)
			Q[i][j] = 0;
	for (int i = 0; i < L; i++)		// заполнение P нулями
		P[i] = 0;

	double **BK = new double*[K*(s - B)*L];	// матрица для оценки bk
	for (int i = 0; i < K*(s - B)*L; i++)
		BK[i] = new double[L];
	for (int i = 0; i < K*(s - B)*L; i++)		// заполнение BK нулями
		for (int j = 0; j < L; j++)
			BK[i][j] = 0;

	char t;
	int l = 0;
	int i = 0;
	int j = 0;
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

	while (true)					// заполнение BK, P частотами
	{
		t = is.get();
		if (is.eof()) break;
		while ((indexof(t) == -1) && (!is.eof()))
			t = is.get();
		if (is.eof()) break;

		i = indexof(t);
		P[i] += 1;

		j = number();
		for (int b = 0; b < s - B; b++)
			BK[j*L*(s - B) + L*b + buffer[b]][i] += 1;
		shiftBuffer(i);
	}

	double maxBlh = -pow(10., 16);
	int bmax = 0;
	double Blh = 0;
	double m = 0;

	for (i = 0; i < K; i++)
	{
		maxBlh = -pow(10., 16);
		bmax = 0;
		for (j = 0; j < s - B; j++)
		{
			Blh = 0;
			for (int l1 = 0; l1 < L; l1++)
			{
				sum = 0;
				for (int l2 = 0; l2 < L; l2++)
					sum += BK[i*L*(s - B) + j*L + l1][l2];
				if (sum > 0)
				{
					for (int l2 = 0; l2 < L; l2++)
						if (BK[i*L*(s - B) + j*L + l1][l2] > 0)
						{
							m = BK[i*L*(s - B) + j*L + l1][l2];
							BK[i*L*(s - B) + j*L + l1][l2] /= sum;
							Blh += m*log(BK[i*L*(s - B) + j*L + l1][l2]);
						}
				}
			}
			if (Blh != 0 && Blh >= maxBlh)
			{
				maxBlh = Blh;
				bmax = j;
			}
		}
		bk[i] = bmax;
		for (int l1 = 0; l1 < L; l1++)
			for (int l2 = 0; l2 < L; l2++)
				Q[i*L + l1][l2] = BK[i*L*(s - B) + bmax*L + l1][l2];
	}

	sum = 0;
	for (i = 0; i < L; i++)
		sum += P[i];
	for (i = 0; i < L; i++)
		P[i] /= sum;

	for (i = 0; i < K*(s - B)*L; i++)
		if (BK[i] != NULL) delete[]BK[i];
	if (BK != NULL) delete[]BK;
}

void ConditionalOrderModel::estimateModelQ(istream &is)
{
	if (is.bad()) cout << "Error! Bad input file for estimating PCModel\n";
	for (int i = 0; i < M*L; i++)		// заполнение Q нулями
		for (int j = 0; j < L; j++)
			Q[i][j] = 0;
	for (int i = 0; i < L; i++)		// заполнение P нулями
		P[i] = 0;

	char t;
	int l = 0;
	int i = 0;
	int j = 0;
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

		j = number();
		Q[mk[j] * L + buffer[bk[j]]][i] += 1;
		shiftBuffer(i);
	}

	for (i = 0; i < M*L; i++)			// нормировка Q, P
	{
		sum = 0;
		for (l = 0; l < L; l++)
			sum += Q[i][l];
		if (sum != 0)
			for (l = 0; l < L; l++)
				Q[i][l] /= sum;
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
}

int ConditionalOrderModel::nextInitialState()
{
	return nextDiscreteRand(P, L);
}

int ConditionalOrderModel::nextState()
{
	int i = number();
	int res = nextDiscreteRand(Q[mk[i] * L + buffer[bk[i]]], L);
	shiftBuffer(res);
	return res;
}

void ConditionalOrderModel::printModel(ostream &os, short full)
{
	MarkovChainModel::printModel(os, full);
	if (full == 0)
		printArray(P, L, os);
	if (full <= 2)
	{
		os << B << endl;
	}
	if (full <= 1)
	{
		os << M << endl;
		if (full != -1)
			printArray(mk, K, os);
		printArray(bk, K, os);
	}
	if (full == 0)
		for (int i = 0; i < M*L; i++)
			printArray(Q[i], L, os);
}
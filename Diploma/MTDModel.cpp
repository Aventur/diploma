#include "MTDModel.h"

double MTDModel::get_buffer_likelihood(Buffer & b)
{
	int h = b[s];
	double res = 0;
	for (int j = 0; j < s; j++)
		res += lambda[j] * Q[b[j]][h];

	return res;
}

double MTDModel::fix_computational_error(double * arr, int size)
{
	int imax = 0;
	double sum = -1.;
	for (int i = 0; i < size; i++)
	{
		sum += arr[i];
		if (arr[i] > arr[imax])
			imax = i;
	}
	arr[imax] -= sum;
	return sum;
}

double MTDModel::vectorIteration(double * target, double * deriv, int size, double eps, double & step, double lh0)
{
	int imin, imax, i;
	double lh, istep;

	imin = 0; imax = 0;
	for (i = 0; i < size; i++)
	{
		if (deriv[i] > deriv[imax])
			imax = i;
		if (deriv[i] < deriv[imin] && target[i] > 0)
			imin = i;
	}

	if (target[imax] < 1 - eps)
	{
		istep = step;
		step = min(target[imin], min(1 - target[imax], step));
		lh = lh0;

		while (step > eps)
		{
			target[imax] += step;
			target[imin] -= step;
			lh = this->likelihood();
			if (lh > lh0)
				break;
			else
			{
				target[imax] -= step;
				target[imin] += step;
				step /= 2;
			}
		}

		if (step == istep || step < eps)
			step *= 2;
		if (lh > lh0)
			return lh;
		else
			return lh0;
	}

	return lh0;
}

double MTDModel::vectorIteration(double * target, double * deriv, int size, double eps, double & step, double lh0, istream * is, int s_pos, int e_pos)
{
	int imin, imax, i;
	double lh, istep;

	imin = 0; imax = 0;
	for (i = 0; i < size; i++)
	{
		if (deriv[i] > deriv[imax])
			imax = i;
		if (deriv[i] < deriv[imin] && target[i] > 0)
			imin = i;
	}

	if (target[imax] < 1 - eps)
	{
		istep = step;
		lh = lh0;
		step = min(target[imin], min(1 - target[imax], step));

		while (step > eps)
		{
			target[imax] += step;
			target[imin] -= step;
			lh = this->likelihood(is, s_pos, e_pos);
			if (lh > lh0)
				break;
			else
			{
				target[imax] -= step;
				target[imin] += step;
				step /= 2;
			}
		}

		if (step == istep || step < eps)
			step *= 2;
		if (lh > lh0)
			return lh;
		else
			return lh0;
	}

	return lh0;
}

void MTDModel::lambda_derivatives(double * dLambda)
{
	Buffer b(L, s + 1);
	ZeroMemory(dLambda, s * sizeof(double));
	int i, j;
	double buffer_lh;
	
	for (i = 0; i < pow(L, s + 1); i++)
	{
		if (nu[i])
		{
			buffer_lh = get_buffer_likelihood(b);
			for (j = 0; j < s; j++)
				dLambda[j] += nu[i] * Q[b[j]][i % L] / buffer_lh;
		}
		b++;
	}
}

void MTDModel::lambda_derivatives(double * dLambda, istream * is, int s_pos, int e_pos)
{
	Buffer b(L, s + 1);
	ZeroMemory(dLambda, s * sizeof(double));
	is->seekg(s_pos, ios::beg);
	int pos, cur, i;
	double buffer_lh;

	for (pos = s_pos; pos < s_pos + s; pos++)
		b.shiftBuffer(indexof(is->get()));

	for (pos = s_pos + s; pos < e_pos; pos++)
	{
		cur = indexof(is->get());
		b.shiftBuffer(cur);
		buffer_lh = get_buffer_likelihood(b);
		for (i = 0; i < s; i++)
			dLambda[i] += Q[b[i]][cur] / buffer_lh;
	}

	refresh_string_stream(is);
}

void MTDModel::Q_derivatives(double * dQl, int l)
{
	Buffer b(L, s + 1);
	ZeroMemory(dQl, L * sizeof(double));
	double sum = 0;
	int i, j, r;

	for (i = 0; i < pow(L, s + 1); i++)
	{
		r = i % L;
		if (!r)
		{
			sum = 0;
			for (j = 0; j < s; j++)
				if (b[j] == l) sum += lambda[j];
		}
		if (sum != 0 && nu[i] != 0)
			dQl[r] += nu[i] * sum / get_buffer_likelihood(b);
		b++;
	}
}

void MTDModel::Q_derivatives(double * dQl, int l, istream * is, int s_pos, int e_pos)
{
	Buffer b(L, s + 1);
	ZeroMemory(dQl, L * sizeof(double));
	is->seekg(s_pos, ios::beg);
	int pos, cur, i;
	double buffer_lh, sum;

	for (pos = s_pos; pos < s_pos + s; pos++)
		b.shiftBuffer(indexof(is->get()));

	for (pos = s_pos + s; pos < e_pos; pos++)
	{
		cur = indexof(is->get());
		b.shiftBuffer(cur);

		sum = 0;
		for (int i = 0; i < s; i++)
			if (b[i] == l)
				sum += lambda[i];
		if (sum)
			dQl[cur] += sum / get_buffer_likelihood(b);
	}

	refresh_string_stream(is);
}

MTDModel::MTDModel()
{
	Q = NULL;
	lambda = NULL;
	P = NULL;
	nu = NULL;
	curDist = NULL;
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
	nu = NULL;
	is.close();
}

MTDModel::MTDModel(int s, int L, char * alpha)
{
	this->s = s;
	this->L = L;
	this->alphabet = new char[L + 1];
	strcpy_s(this->alphabet, L + 1, alpha);
	this->alphabet[L] = '\0';
	//this->buffer = new int[L];

	lambda = new double[s];
	for (int i = 0; i < s; i++)
		lambda[i] = 1. / s;
	Q = new double*[L];
	for (int i = 0; i < L; i++)
	{
		Q[i] = new double[L];
		for (int j = 0; j < L; j++)
			Q[i][j] = 1. / L;
	}

	P = NULL;
	nu = NULL;
	curDist = NULL;
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
	if (curDist != NULL) delete[]curDist;
	if (nu != NULL) delete[]nu;
}

void MTDModel::estimateNu(istream * is)
{
	is->seekg(0, ios::end);
	int T = is->tellg();
	is->seekg(0, ios::beg);
	estimateNu(is, 0, T);
}

void MTDModel::estimateNu(istream * is, int s_pos, int e_pos)
{
	if (nu == NULL) nu = new unsigned int[pow(L, s + 1)];
	ZeroMemory(nu, sizeof(unsigned int) * pow(L, s + 1));
	is->seekg(s_pos, ios::beg);

	Buffer *b = new Buffer(L, s + 1);
	int pos;
	for (pos = s_pos; pos < s_pos + s; pos++)
		b->shiftBuffer(indexof(is->get()));

	for (pos = s_pos + s; pos < e_pos; pos++)
	{
		b->shiftBuffer(indexof(is->get()));
		nu[b->getNumber()]++;
	} 

	delete b;
}

void MTDModel::estimateInitialParameters(istream * is)
{
	is->seekg(0, ios::end);
	int T = is->tellg();
	is->seekg(0, ios::beg);
	estimateInitialParameters(is, 0, T);
}

void MTDModel::estimateInitialParameters(istream * is, int s_pos, int e_pos)
{
	is->seekg(s_pos, ios::beg);
	int T = e_pos - s_pos;

	if (P == NULL) P = new double[L];
	ZeroMemory(P, sizeof(double) * L);

	int i, j, k, t, cur;		// выделение памяти на вспомогательные переменные
	double ***Pi = new double**[s];
	double ***Z = new double**[s];
	for (i = 0; i < s; i++)
	{
		Pi[i] = new double*[L];
		Z[i] = new double*[L];
		for (j = 0; j < L; j++)
		{
			Pi[i][j] = new double[L];
			Z[i][j] = new double[L];
			ZeroMemory(Pi[i][j], sizeof(double) * L);
		}
	}
	double **D = new double*[L];
	for (i = 0; i < L; i++)
		D[i] = new double[L];

	Buffer b(L, s + 1);			// подсчёт абсолютных частот
	for (t = 1; t <= T; t++)
	{
		cur = indexof(is->get());

		b.shiftBuffer(cur);
		if (t >= s + 1 && t <= T - s + 1)
			P[cur] += 1.;
		for (j = 1; j <= s; j++)
			if (t >= s + j && t <= T - s + j)
				Pi[j - 1][b[j - 1]][cur] += 1.;
	}

	for (i = 0; i < L; i++)		// подсчёт относительных частот
		P[i] /= T - 2 * s + 1;
	for (j = 0; j < s; j++)
		for (i = 0; i < L; i++)
			for (k = 0; k < L; k++)
				Pi[j][i][k] /= T - 2 * s + 1;

	// вычисление оценки Q
	for (i = 0; i < L; i++)
		ZeroMemory(Q[i], sizeof(double) * L);
	for (j = 0; j < s; j++)
		for (k = 0; k < L; k++)
			for (i = 0; i < L; i++)
				Q[k][i] += Pi[j][k][i];
	for (k = 0; k < L; k++)
		for (i = 0; i < L; i++)
			Q[k][i] = (P[k] > 0) ? (Q[k][i] / P[k] - (s - 1) * P[i]) : (1. / L);

	for (j = 0; j < s; j++)		// заполнение Z
		for (k = 0; k < L; k++)
			for (i = 0; i < L; i++)
				Z[j][k][i] = (P[k] > 0) ? (Pi[j][k][i] / P[k] - P[i]) : (0.);

	for (k = 0; k < L; k++)		// заполнение D
		for (i = 0; i < L; i++)
			D[k][i] = Q[k][i] - P[i];

	// вычисление оценки lambda
	ZeroMemory(lambda, sizeof(double) * s);
	double sumD2 = 0;
	for (k = 0; k < L; k++)
		for (i = 0; i < L; i++)
			sumD2 += D[k][i] * D[k][i];
	for (j = 0; j < s; j++)
	{
		for (k = 0; k < L; k++)
			for (i = 0; i < L; i++)
				lambda[j] += Z[j][k][i] * D[k][i];
		lambda[j] /= sumD2;
	}

	// исправление погрешности (для выполнения условия нормировки)
	fix_computational_error(lambda, s);
	for (i = 0; i < L; i++)
		fix_computational_error(Q[i], L);

	// коррекция результата (проекция оценок на симплекс)
	project_to_simplex(lambda, s);
	for (i = 0; i < L; i++)
		project_to_simplex(Q[i], L);

	for (j = 0; j < s; j++)		// освобождение памяти
	{
		for (i = 0; i < L; i++)
		{
			if (Pi[j][i] != NULL) delete[]Pi[j][i];
			if (Z[j][i] != NULL) delete[]Z[j][i];
		}
		if (Pi[j] != NULL) delete[]Pi[j];
		if (Z[j] != NULL) delete[]Z[j];
	}
	if (Pi != NULL) delete[]Pi;
	if (Z != NULL) delete[]Z;
	for (i = 0; i < L; i++)
		if (D[i] != NULL) delete[]D[i];
	if (D != NULL) delete[]D;
	refresh_string_stream(is);
}

void MTDModel::printModel(ostream &os)
{
	MarkovChainModel::printModel(os);
	printArray(P, L, os);
	printArray(lambda, s, os);
	for (int i = 0; i < L; i++)
		printArray(Q[i], L, os);
}

double MTDModel::iterativeEstimation(ostream & os, double eps, int max_iter, double start_step, bool report)
{
	double lh0, lh, lstep;
	double *qsteps = new double[L];
	int i, iter;

	double **dQ = new double*[L];
	for (i = 0; i < L; i++)
		dQ[i] = new double[L];
	double *dLambda = new double[s];

	lstep = start_step;
	for (i = 0; i < L; i++)
		qsteps[i] = start_step;
	lh0 = -DBL_MAX;
	lh = likelihood();

	iter = 0;
	while (lh > lh0 && i <= max_iter)
	{
		if (report)
		{
			os << "\n\nIteration " << ++iter << ":\n";
			cout << endl << "Iteration " << iter << ":\n";
			printModel(cout);
		}

		lh0 = lh;
		lambda_derivatives(dLambda);

		if (report)
		{
			os << "Lambda iteration:\n";
			printArray(dLambda, L, os);
			os << fixed << "lh before: " << lh0 << ", step before: " << lstep
				<< endl;
		}

		lh = vectorIteration(lambda, dLambda, s, eps, lstep, lh0);

		if (report)
			os << fixed << "lh after: " << lh << ", step after: " << lstep
			<< endl;

		for (i = 0; i < L; i++)
		{
			Q_derivatives(dQ[i], i);
			if (report)
			{
				os << "Q" << i << " iteration:\n";
				printArray(dQ[i], L, os);
				os << fixed << "lh before: " << lh << ", step before: " << qsteps[i]
					<< endl;
			}

			lh = vectorIteration(Q[i], dQ[i], L, eps, qsteps[i], lh);

			if (report)
				os << fixed << "lh after: " << lh << ", step after: " << qsteps[i]
				<< endl;
		}

		if (report)
		{
			printModel(os);
			os << endl << "-------------------------------" << endl;
		}
	}

	if (qsteps != NULL) delete[]qsteps;
	for (i = 0; i < L; i++)
		if (dQ[i] != NULL) delete[]dQ[i];
	if (dQ != NULL) delete[]dQ;
	if (dLambda != NULL) delete[]dLambda;
	return lh;
}

double MTDModel::iterativeEstimation_Nu(istream * is, int s_pos, int e_pos, ostream &os, double eps, int max_iter, double start_step, bool report)
{
	estimateNu(is, s_pos, e_pos);
	return iterativeEstimation(os, eps, max_iter, start_step, report);
}

double MTDModel::iterativeEstimation(istream * is, int s_pos, int e_pos, ostream & os,
	double eps, int max_iter, double start_step, bool report)
{
	double lh0, lh, lstep;
	double *qsteps = new double[L];
	int i, iter;

	double **dQ = new double*[L];
	for (i = 0; i < L; i++)
		dQ[i] = new double[L];
	double *dLambda = new double[s];

	lstep = start_step;
	for (i = 0; i < L; i++)
		qsteps[i] = start_step;
	lh0 = -DBL_MAX;
	lh = likelihood(is, s_pos, e_pos);

	iter = 0;
	while (lh > lh0 && i <= max_iter)
	{
		if (report)
		{
			os << "\n\nIteration " << ++iter << ":\n";
			cout << endl << "Iteration " << iter << ":\n";
			printModel(cout);
		}

		lh0 = lh;
		lambda_derivatives(dLambda, is, s_pos, e_pos);

		if (report)
		{
			os << "Lambda iteration:\n";
			printArray(dLambda, L, os);
			os << fixed << "lh before: " << lh0 << ", step before: " << lstep
				<< endl;
		}

		lh = vectorIteration(lambda, dLambda, s, eps, lstep, lh0, is, s_pos, e_pos);

		if (report)
			os << fixed << "lh after: " << lh << ", step after: " << lstep
				<< endl;

		for (i = 0; i < L; i++)
		{
			Q_derivatives(dQ[i], i, is, s_pos, e_pos);
			if (report)
			{
				os << "Q" << i << " iteration:\n";
				printArray(dQ[i], L, os);
				os << fixed << "lh before: " << lh << ", step before: " << qsteps[i]
					<< endl;
			}

			lh = vectorIteration(Q[i], dQ[i], L, eps, qsteps[i], lh, is, s_pos, e_pos);

			if (report)
				os << fixed << "lh after: " << lh << ", step after: " << qsteps[i]
				<< endl;
		}

		if (report)
		{
			printModel(os);
			os << endl << "-------------------------------" << endl;
		}
	}

	if (qsteps != NULL) delete[]qsteps;
	for (i = 0; i < L; i++)
		if (dQ[i] != NULL) delete[]dQ[i];
	if (dQ != NULL) delete[]dQ;
	if (dLambda != NULL) delete[]dLambda;
	return lh;
}

double MTDModel::likelihood()
{
	double lh = 0;
	Buffer b(L, s + 1);
	int i;

	for (i = 0; i < pow(L, s + 1); i++)
	{
		if (nu[i])
			lh += nu[i] * log(get_buffer_likelihood(b));
		b++;
	}

	return lh;
}

double MTDModel::likelihood(istream * is)
{
	double lh = 0;
	int cur;
	Buffer b(L, s+1);

	for (int i = 0; i < s; i++)
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

double MTDModel::likelihood(istream * is, int s_pos, int e_pos)
{
	double lh = 0;
	int cur, pos;
	Buffer b(L, s + 1);
	is->seekg(s_pos, ios::beg);

	for (pos = s_pos; pos < s_pos + s; pos++)
	{
		cur = indexof(is->get());
		b.shiftBuffer(cur);
	}

	for (pos = s_pos + s; pos < e_pos; pos++)
	{
		cur = indexof(is->get());
		if (cur == -1)
			continue;
		b.shiftBuffer(cur);
		double buffer_lh = get_buffer_likelihood(b);	/// !check
		if (buffer_lh)
			lh += log(buffer_lh);
	}

	refresh_string_stream(is);
	return lh;
}

__int64 MTDModel::numberOfParams()
{
	return L * (L - 1) + s - 1;
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
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
	os << endl;
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
				Q[i - j] = sum?(nu[i - j] / sum): 1./L;
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
	ZeroMemory(C, sizeof(double) * s);

	double c = 0;					// вычисление C
	for (int i = 0; i < L; i++)
		c += P[i] * P[i];
	c *= 2 * pow(L, s - 1) * ro * ro;

	do
	{
		for (int j = 0; j < s; j++)
			C[j] += P[b[j]] - Q[L*b.getNumber() + b[j]];
		b++;
	} while (b.getNumber());
	for (int i = 0; i < s; i++)
		C[i] = c - 2 * ro * C[i];

	double prod;
	for (int i = 0; i < s; i++)		// вычисление lambda
	{
		prod = 0;
		for (int j = 0; j < s; j++)
			prod += (s * ((s-i-1) == j) - 1) * C[j];
		lambda[i] = prod / (2 * s * (L - 1) * pow(L, s - 1) * ro * ro) + 1. / s;
	}

	project_to_simplex(lambda, s);
	return lambda;
}

double JacobsLewisModel::get_buffer_likelihood(Buffer &b)
{
	int h = b[s];
	double res = 0;
	for (int j = 0; j < s; j++)
		if (h == b[j])
			res += lambda[s - j - 1];

	return ((1. - ro) * P[h] + ro * res);
}

void JacobsLewisModel::P_derivatives(double * d_P)
{
	int i;
	ZeroMemory(d_P, sizeof(double) * L);

	Buffer b(L, s + 1);
	for (i = 0; i < pow(L, s + 1); i++)
	{
		if (nu[i])
			d_P[b[s]] += nu[i] / get_buffer_likelihood(b);
		b++;
	}

	for (i = 0; i < L; i++)
		d_P[i] *= (1. - ro);
}

void JacobsLewisModel::Lambda_derivatives(double * d_lambda)
{	// with use of nu
	int i, k, h;
	ZeroMemory(d_lambda, sizeof(double) * s);

	Buffer b(L, s + 1);
	double blh = 0;
	for (i = 0; i < pow(L, s + 1); i++)
	{
		if (nu[i])
		{
			blh = 0;
			h = b[s];
			for (k = 0; k < s; k++)
				if (h == b[k])
				{
					if (!blh) blh = get_buffer_likelihood(b);
					d_lambda[s - k - 1] += nu[i] / blh;
				}
		}
		b++;
	}

	for (i = 0; i < s; i++)
		d_lambda[i] *= ro;
}

double JacobsLewisModel::Ro_derivative()
{
	int i, j, h;
	double d_ro = 0;

	Buffer b(L, s + 1);
	double sum = 0;
	for (i = 0; i < pow(L, s + 1); i++)
	{
		if (nu[i])
		{
			h = b[s];
			sum = 0;
			for (j = 0; j < s; j++)
				if (h == b[j])
					sum += lambda[s - j - 1];

			d_ro += nu[i] * (sum - P[h]) / ((1. - ro) * P[h] + ro * sum);
		}
		b++;
	}

	return d_ro;
}

void JacobsLewisModel::all_derivatives(double * d_P, double * d_lambda, double & d_ro)
{
	int i, j, h;
	ZeroMemory(d_P, sizeof(double) * L);
	ZeroMemory(d_lambda, sizeof(double) * s);
	d_ro = 0;

	Buffer b(L, s + 1);
	double sum, blh;
	for (i = 0; i < pow(L, s + 1); i++)
	{
		if (nu[i])
		{
			h = b[s];
			sum = 0;
			for (j = 0; j < s; j++)
				if (h == b[j])
					sum += lambda[s - j - 1];
			blh = ((1. - ro) * P[h] + ro * sum);

			for (j = 0; j < s; j++)
				if (h == b[j])
					d_lambda[s - j - 1] += nu[i] / blh;
			d_P[h] += nu[i] / blh;
			d_ro += nu[i] * (sum - P[h]) / blh;
		}
		b++;
	}
	for (i = 0; i < L; i++)
		d_P[i] *= (1. - ro);
	for (i = 0; i < s; i++)
		d_lambda[i] *= ro;
}

void JacobsLewisModel::P_derivatives(double * d_P, istream * is, int s_pos, int e_pos)
{
	int i, c;
	ZeroMemory(d_P, sizeof(double) * L);

	Buffer b(L, s + 1);
	is->seekg(s_pos, ios::beg);
	for (i = 0; i < s; i++)			// а надо ли?
	{
		c = indexof(is->get());
		d_P[c] += 1.;
		b.shiftBuffer(c);
	}

	for (i = s_pos + s; i < e_pos; i++)
	{
		c = indexof(is->get());
		b.shiftBuffer(c);
		d_P[c] += 1. / get_buffer_likelihood(b);
	}

	for (i = 0; i < L; i++)
		d_P[i] *= (1. - ro);
	refresh_string_stream(is);
}

void JacobsLewisModel::Lambda_derivatives(double * d_lambda, istream * is, int s_pos, int e_pos)
{
	int i, c, k;
	ZeroMemory(d_lambda, sizeof(double) * s);

	Buffer b(L, s + 1);
	is->seekg(s_pos, ios::beg);
	for (i = 0; i < s; i++)	
		b.shiftBuffer(indexof(is->get()));

	double blh;
	for (i = s_pos + s; i < e_pos; i++)
	{
		b.shiftBuffer(indexof(is->get()));

		blh = 0;
		c = b[s];
		for (k = 0; k < s; k++)
			if (c == b[k])
			{
				if (!blh) blh = get_buffer_likelihood(b);
				d_lambda[s - k - 1] += 1. / blh;
			}
	}

	for (i = 0; i < s; i++)
		d_lambda[i] *= ro;
	refresh_string_stream(is);
}

double JacobsLewisModel::Ro_derivative(istream * is, int s_pos, int e_pos)
{
	int i, j, h;
	double d_ro = 0;

	Buffer b(L, s + 1);
	is->seekg(s_pos, ios::beg);
	for (i = 0; i < s; i++)
		b.shiftBuffer(indexof(is->get()));

	double sum = 0;
	for (i = s_pos + s; i < e_pos; i++)
	{
		h = indexof(is->get());
		b.shiftBuffer(h);
		sum = 0;
		for (j = 0; j < s; j++)
			if (h == b[j])
				sum += lambda[s - j - 1];

		d_ro += (sum - P[h]) / ((1. - ro) * P[h] + ro * sum);
	}

	refresh_string_stream(is);
	return d_ro;
}

void JacobsLewisModel::all_derivatives(double * d_P, double * d_lambda, double & d_ro, istream * is, int s_pos, int e_pos)
{
	int i, j, h;
	ZeroMemory(d_P, sizeof(double) * L);
	ZeroMemory(d_lambda, sizeof(double) * s);
	d_ro = 0;

	Buffer b(L, s + 1);
	is->seekg(s_pos, ios::beg);
	for (i = 0; i < s; i++)	
	{
		h = indexof(is->get());
		d_P[h] += 1.;
		b.shiftBuffer(h);
	}

	double sum, blh;
	for (i = s_pos + s; i < e_pos; i++)
	{
		h = indexof(is->get());
		b.shiftBuffer(h);
		sum = 0;
		for (j = 0; j < s; j++)
			if (h == b[j])
				sum += lambda[s - j - 1];
		blh = ((1. - ro) * P[h] + ro * sum);

		for (j = 0; j < s; j++)
			if (h == b[j])
				d_lambda[s - j - 1] += 1. / blh;
		d_P[h] += 1. / blh;
		d_ro += (sum - P[h]) / blh;
	}
	for (i = 0; i < L; i++)
		d_P[i] *= (1. - ro);
	for (i = 0; i < s; i++)
		d_lambda[i] *= ro;
	refresh_string_stream(is);
}

double JacobsLewisModel::vectorIteration(double * target, double * deriv, int size, 
	double eps, double & step, double lh0)
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

double JacobsLewisModel::singleIteration(double & target, double deriv, 
	double eps, double & step, double lh0)
{
	if (abs(deriv) < eps)
		return lh0;
	double lh, istep;
	double direction = (deriv > 0) ? (1) : (-1);

	istep = step;
	step = min(step, (deriv > 0)?(1 - target):(target));

	while (step > eps)
	{
		target += step * direction;
		lh = this->likelihood();
		if (lh > lh0)
			break;
		else
		{
			target -= step * direction;
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

void JacobsLewisModel::iterativeEstimation(istream * is, ostream &os, double eps)
{
	estimateInitialParameters(is);//estimateNu(is);
	refresh_string_stream(is);
	double lh0, lh, lstep, pstep, rstep;
	int i = 0;

	double *dP = new double[L];
	double *dLambda = new double[s];
	double dRo = 0;

	lstep = 0.1; pstep = 0.1; rstep = 0.1;
	lh0 = -DBL_MAX;
	lh = likelihood(is);

	while (lh > lh0)
	{
		os << endl << "Iteration " << ++i << ":\n";
		cout << endl << "Iteration " << i << ":\n";
		printModel(cout);

		lh0 = lh;
		P_derivatives(dP);

		os << "P iteration:\n";
		printModel(os);
		printArray(dP, L, os);
		os << fixed << lh0 << " " << pstep << endl << "--------" << endl;

		lh = vectorIteration(P, dP, L, eps, pstep, lh0);

		printModel(os);
		os << fixed << lh << " " << pstep << endl << "--------" << endl;

		Lambda_derivatives(dLambda);

		os << "Lambda iteration:\n";
		printArray(dLambda, s, os);
		os << fixed << lh << " " << lstep << endl << "--------" << endl;

		lh = vectorIteration(lambda, dLambda, s, eps, lstep, lh);

		printModel(os);
		os << fixed << lh << " " << lstep << endl << "--------" << endl;

		dRo = Ro_derivative();

		os << "Ro iteration:\n";
		os << dRo << endl;
		os << fixed << lh << " " << rstep << endl << "--------" << endl;

		lh = singleIteration(ro, dRo, eps, rstep, lh);

		printModel(os);
		os << fixed << lh << " " << rstep << endl << "--------" << endl;
	}

	if (dP != NULL) delete[]dP;
	if (dLambda != NULL) delete[]dLambda;
}

double JacobsLewisModel::vectorIteration(double * target, double * deriv, int size, double eps, double & step, double lh0, istream * is, int s_pos, int e_pos)
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

double JacobsLewisModel::singleIteration(double & target, double deriv, double eps, double & step, double lh0, istream * is, int s_pos, int e_pos)
{
	if (abs(deriv) < eps)
		return lh0;
	double lh, istep;
	lh = lh0;
	double direction = (deriv > 0) ? (1) : (-1);

	istep = step;
	step = min(step, (deriv > 0) ? (1 - target) : (target));

	while (step > eps)
	{
		target += step * direction;
		lh = this->likelihood(is, s_pos, e_pos);
		if (lh > lh0)
			break;
		else
		{
			target -= step * direction;
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

double JacobsLewisModel::iterativeEstimation(istream * is, int s_pos, int e_pos, double eps, int max_iter)
{
	double lh0, lh, lstep, pstep, rstep;
	int i = 0;

	double *dP = new double[L];
	double *dLambda = new double[s];
	double dRo = 0;

	lstep = 0.1; pstep = 0.1; rstep = 0.1;
	lh0 = -DBL_MAX;
	lh = likelihood(is, s_pos, e_pos);

	while (lh > lh0 && i < max_iter)
	{
		i++;
		lh0 = lh;

		//all_derivatives(dP, dLambda, dRo, is, s_pos, e_pos);
		P_derivatives(dP, is, s_pos, e_pos);
		lh = vectorIteration(P, dP, L, eps, pstep, lh0, is, s_pos, e_pos);

		Lambda_derivatives(dLambda, is, s_pos, e_pos);
		lh = vectorIteration(lambda, dLambda, s, eps, lstep, lh, is, s_pos, e_pos);

		dRo = Ro_derivative(is, s_pos, e_pos);
		lh = singleIteration(ro, dRo, eps, rstep, lh, is, s_pos, e_pos);
	}

	if (dP != NULL) delete[]dP;
	if (dLambda != NULL) delete[]dLambda;
	return lh;
}

void JacobsLewisModel::iterativeEstimation(istream * is, int s_pos, int e_pos, ostream & os, double eps)
{
	double lh0, lh, lstep, pstep, rstep;
	int i = 0;

	double *dP = new double[L];
	double *dLambda = new double[s];
	double dRo = 0;

	lstep = 0.1; pstep = 0.1; rstep = 0.1;
	lh0 = -DBL_MAX;
	lh = likelihood(is);

	while (lh > lh0)
	{
		os << endl << "Iteration " << ++i << ":\n";
		cout << endl << "Iteration " << i << ":\n";
		printModel(cout);

		lh0 = lh;
		P_derivatives(dP, is,  s_pos, e_pos);
		//all_derivatives(dP, dLambda, dRo, is, s_pos, e_pos);

		os << "P iteration:\n";
		printModel(os);
		printArray(dP, L, os);
		os << fixed << lh0 << " " << pstep << endl << "--------" << endl;

		lh = vectorIteration(P, dP, L, eps, pstep, lh0, is, s_pos, e_pos);

		printModel(os);
		os << fixed << lh << " " << pstep << endl << "--------" << endl;

		Lambda_derivatives(dLambda, is, s_pos, e_pos);

		os << "Lambda iteration:\n";
		printArray(dLambda, s, os);
		os << fixed << lh << " " << lstep << endl << "--------" << endl;

		lh = vectorIteration(lambda, dLambda, s, eps, lstep, lh, is, s_pos, e_pos);

		printModel(os);
		os << fixed << lh << " " << lstep << endl << "--------" << endl;

		dRo = Ro_derivative(is, s_pos, e_pos);

		os << "Ro iteration:\n";
		os << dRo << endl;
		os << fixed << lh << " " << rstep << endl << "--------" << endl;

		lh = singleIteration(ro, dRo, eps, rstep, lh, is, s_pos, e_pos);

		printModel(os);
		os << fixed << lh << " " << rstep << endl << "--------" << endl;
	}

	if (dP != NULL) delete[]dP;
	if (dLambda != NULL) delete[]dLambda;
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

	refresh_string_stream(is);
	//if (frequencies != NULL) delete[]frequencies;
}

void JacobsLewisModel::estimateInitialParameters(istream * is)
{
	estimateNu(is);
	refresh_string_stream(is);

	// вычислим оценку матрицы Q
	double *Q = new double[pow(L, s + 1)];
	estimateQ(Q);

	estimateP(is);				// оценка начального значения P
	refresh_string_stream(is);

	estimateInitialRo(Q);		// оценка начального значения Ro
	estimateInitialLambda(Q);	// оценка начального значения Lambda

	if (Q != NULL) delete[]Q;
}

void JacobsLewisModel::setRandomInitialParameters()
{
	generateRandomDistribution(lambda, s);
	generateRandomDistribution(P, L);
	ro = bsv->next();
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

double JacobsLewisModel::likelihood(istream *is)
{
	double lh = 0;
	int cur;
	Buffer b(L, s);

	for (int i = 0; i < s; i++)
	{
		cur = indexof(is->get());
		b.shiftBuffer(cur);
		lh += log(P[cur]);
	}

	double sum;
	int i;
	for (char c = is->get(); !is->eof(); c = is->get())
	{
		cur = indexof(c);
		if (cur == -1)
			continue;

		sum = 0;
		for (i = 0; i < s; i++)
			if (cur == b[i])
				sum += lambda[s - i - 1];
		lh += log((1 - ro) * P[cur] + ro * sum);

		b.shiftBuffer(cur);
	}

	refresh_string_stream(is);
	return lh;
}

double JacobsLewisModel::likelihood(istream *is, int s_pos, int e_pos)
{
	double lh = 0;
	int cur;
	Buffer b(L, s);
	is->seekg(s_pos, ios::beg);
	int pos;

	for (pos = s_pos; pos < s_pos + s; pos++)
	{
		cur = indexof(is->get());
		b.shiftBuffer(cur);
		lh += log(P[cur]);
	}

	double sum;
	int i;
	for (char c = is->get(); pos < e_pos; c = is->get())
	{
		cur = indexof(c);
		if (cur == -1)
			continue;

		sum = 0;
		for (i = 0; i < s; i++)
			if (cur == b[i])
				sum += lambda[s - i - 1];
		lh += log((1 - ro) * P[cur] + ro * sum);

		b.shiftBuffer(cur);
		pos++;
	}

	refresh_string_stream(is);
	return lh;
}

__int64 JacobsLewisModel::numberOfParams()
{
	return (L + s - 1);
}
#include "MCModel.h"

#ifndef JLMODEL
#define JLMODEL
/*
Формат входного файла:
порядок модели
число состояний
алфавит
распределение пи через пробел
распределение лямбда через пробел
число ро
*/
class JacobsLewisModel : public MarkovChainModel
{
	double ro;			// значение ro
	double *lambda;		// вектор вероятностей lambda
	double *P;			// вектор вероятностей pi

	// опциональные поля (для жадного вычисления оценок)
	unsigned int *nu;			// массив для хранения матрицы частот
	void estimateNu(istream *is); 
	double *estimateQ(double *Q); // после оценки Nu

	double estimateInitialRo(double *Q); // после оценки P
	double *estimateInitialLambda(double *Q); // после после P

	double get_buffer_likelihood(Buffer b);// вспомогательная для следующих 6
	void P_derivatives(double* d_P);	// после оценки nu, P, lambda, ro
	void Lambda_derivatives(double* d_lambda);
	double Ro_derivative();
	void all_derivatives(double *d_P, double *d_lambda, double &d_ro); 

	void P_derivatives(double* d_P, istream *is, int s_pos, int e_pos);
	void Lambda_derivatives(double* d_lambda, istream *is, int s_pos, int e_pos);
	double Ro_derivative(istream *is, int s_pos, int e_pos);
	void all_derivatives(double *d_P, double *d_lambda, double &d_ro, 
		istream *is, int s_pos, int e_pos); // если нужно посчитать всё за один проход

	__int64 next1(Buffer &b, int i); // s+1-я компонента должна быть равна i, остальные != i
	__int64 next2(Buffer &b, int i); // все компоненты равны i

public:
	JacobsLewisModel();
	JacobsLewisModel(int s, int L, char* alpha);	// пустая модель, для оценки параметров;
													// lambda = (1/s), P = (1/L), ro = 0.5
	JacobsLewisModel(char *fname);					// инициализировать модель параметрами из файла
	~JacobsLewisModel();

	void estimateP(istream *is);
	void estimateInitialParameters(istream *is);	// !после оценки nu и P

	double likelihood(istream *is);
	__int64 numberOfParams();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
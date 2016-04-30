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

	double likelihood();
	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
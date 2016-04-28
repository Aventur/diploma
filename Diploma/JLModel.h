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

	unsigned int *nu;			// массив для хранения матрицы частот

public:
	JacobsLewisModel();
	JacobsLewisModel(int s, int L, char* alpha);	// пустая модель, для оценки параметров;
													// lambda = (1/s), P = (1/L), ro = 0.5
	JacobsLewisModel(char *fname);					// инициализировать модель параметрами из файла
	~JacobsLewisModel();

	void estimateNu(istream &is);
	void estimateP(istream &is);

	double likelihood();
	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
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

public:
	JacobsLewisModel();
	JacobsLewisModel(char *fname);
	~JacobsLewisModel();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
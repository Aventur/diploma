#include "MCModel.h"

#ifndef MTDMODEL
#define MTDMODEL
/*
Формат входного файла:
порядок модели
число состояний
алфавит
распределение пи через пробел
веса лямбда через пробел
матрица Q по строкам
*/
class MTDModel : public MarkovChainModel
{
	double **Q;			// матрица Q
	double *lambda;		// вектор весов лямбда
	double *P;			// вектор одномерного стационарного распределения

	double *curDist;	// вспомогательный вектор

public:
	MTDModel();
	MTDModel(char *fname);
	~MTDModel();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
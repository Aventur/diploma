#include "MCModel.h"

#ifndef PCMODEL
#define PCMODEL
/*
Формат входного файла: (соответствующее значение full)
порядок модели (0, 1, 2)
число состояний (all)
алфавит (all)
распределение пи через пробел (0)
число связей r (0, 1, 2)
шаблон связей (r чисел из 0..s-1) (0, 1)
матрица Q по строкам из L чисел (0)
*/
class PartialConnectionsModel : public MarkovChainModel
{
	double *P;			// одномерное стационарное распределение
	double **Q;			// матрица переходов в двумерной форме
	int r;				// число связей
	int *m;				// шаблон связей
	int c;				// число строк в матрице Q

	int number();		// возвращает номер текущей подпоследовательности в буфере с учётом шаблона
	bool nextM(int level);

public:
	PartialConnectionsModel();
	PartialConnectionsModel(char *fname, short full = 0);
	~PartialConnectionsModel();

	void setParam(int p1, int p2);
	void estimateModel(istream &is);	// известны только s, r, L, alpha (full = 2)
	double estimateModelQ(istream &is);	// известно всё, кроме Q (full = 1)
	int nextInitialState();
	int nextState();
	void printModel(ostream &os, short full);
};

#endif
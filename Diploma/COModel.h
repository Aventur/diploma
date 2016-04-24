#include "MCModel.h"

#ifndef COMODEL
#define COMODEL
/*
Формат входного файла:
порядок модели (0, 1, 2)
число состояний (all)
алфавит (all)
распределение пи через пробел (0)
длина БФП (одно число) (0, 1, 2)
число M различных матриц Q (одно число) (0, 1)
номера матриц, соответствующих БФП (K чисел mk) (0, 1)
номера состояний, из которых осуществляется переход, для каждого БФП (K чисел bk) (0, 1)
M матриц Q по строкам (0)
*/
class ConditionalOrderModel : public MarkovChainModel
{
	double *P;			// стационарное распределение
	int B;				// длина БФП
	int M;				// число различных матриц Q
	int K;				// число различных вариантов БФП
	int *mk;			// номера матриц, соответствующих БФП (начиная с 0)
	int *bk;			// номера состояний, из которых осуществляется переход, соотв. БФП (начиная с 0)
	double **Q;			// M матриц в одной, по строкам

	int number();		// номер текущего БФП (того, что в буфере)

public:
	ConditionalOrderModel();
	ConditionalOrderModel(char *fname, short full = 0);
	~ConditionalOrderModel();

	void setParam(int p1, int p2);
	void estimateModel(istream &is);
	void estimateModelB(istream &is);
	void estimateModelQ(istream &is);	// для заданных alphabet, L, s, B, M, mk, bk
	int nextInitialState();
	int nextState();
	void printModel(ostream &os, short full = 0);
};

#endif
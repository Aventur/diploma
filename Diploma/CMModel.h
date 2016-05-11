#include "MCModel.h"

#ifndef CLASSICMARKOVMODEL
#define CLASSICMARKOVMODEL
/*
‘ормат входного файла:
пор€док модели
число состо€ний
алфавит
матрица P по строкам
*/
class ClassicMarkovModel : public MarkovChainModel
{
public:
	double *Q;					// матрица веро€тностей переходов, записанна€ в строку
	unsigned int *nu;			// массив дл€ хранени€ матрицы частот

	void estimateNu(istream *is);
	void estimateNu(istream *is, int s_pos, int e_pos);

	double get_buffer_likelihood(Buffer &b);// вспомогательна€

//public:
	ClassicMarkovModel();
	ClassicMarkovModel(int s, int L, char* alpha);	// пуста€ модель, дл€ оценки параметров;
	ClassicMarkovModel(char *fname);	// инициализировать модель параметрами из файла
	~ClassicMarkovModel();

	double estimateQ(istream *is);
	/*returns likelihood of the estimated model*/
	double estimateQ(istream *is, int s_pos, int e_pos);

	double likelihood(istream *is);
	double likelihood(istream *is, int s_pos, int e_pos);
	__int64 numberOfParams();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
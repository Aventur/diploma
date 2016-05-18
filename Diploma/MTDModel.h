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

	double *curDist;	// вспомогательный вектор для генерации

	unsigned int *nu;	// матрица частот для быстрого вычисления правдоподобия
						// и производных по длинным файлам

	double get_buffer_likelihood(Buffer &b);
	/* обеспечение условия нормировки для дискретного распределения */
	double fix_computational_error(double *arr, int size);

	/* при заданной nu, для оч. больших T и малых s*/
	double vectorIteration(double *target, double *deriv, int size,
		double eps, double &step, double lh0);
	/* likelihood вычисляется по последовательности */
	double vectorIteration(double *target, double *deriv, int size,
		double eps, double &step, double lh0, istream *is, int s_pos, int e_pos);

	/*derivatives !!! по две функции для лямбда, строки Q 
	(и всех сразу - пока не надо)*/
	void lambda_derivatives(double *dLambda);	// по nu
	void lambda_derivatives(double *dLambda, istream *is, int s_pos, int e_pos);

	/* производные строки матрицы Q с номером l*/
	void Q_derivatives(double *dQl, int l);		// по nu
	void Q_derivatives(double *dQl, int l, istream *is, int s_pos, int e_pos);

	/* с оценкой nu */
	double iterativeEstimation(ostream &os, double eps = 0.001, int max_iter = 100,
		double start_step = 0.1, bool report = 0);

public:
	MTDModel();
	MTDModel(char *fname);					// для генерации
	MTDModel(int s, int L, char *alpha);	// для оценки параметров, с установкой по умолчанию
	~MTDModel();

	void estimateNu(istream *is);
	void estimateNu(istream *is, int s_pos, int e_pos);

	void estimateInitialParameters(istream *is);
	void estimateInitialParameters(istream *is, int s_pos, int e_pos);

	/* с оценкой nu */
	double iterativeEstimation_Nu(istream *is, int s_pos, int e_pos, ostream &os,
		double eps = 0.001, int max_iter = 100,	double start_step = 0.1, bool report = 0);
	/* без оценки nu */
	double iterativeEstimation(istream *is, int s_pos, int e_pos, ostream &os,
		double eps = 0.001, int max_iter = 100, double start_step = 0.1, bool report = 0);

	double likelihood();
	double likelihood(istream *is);
	double likelihood(istream *is, int s_pos, int e_pos);
	__int64 numberOfParams();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
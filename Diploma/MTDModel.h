#include "MCModel.h"

#ifndef MTDMODEL
#define MTDMODEL
/*
������ �������� �����:
������� ������
����� ���������
�������
������������� �� ����� ������
���� ������ ����� ������
������� Q �� �������
*/
class MTDModel : public MarkovChainModel
{
	double **Q;			// ������� Q
	double *lambda;		// ������ ����� ������
	double *P;			// ������ ����������� ������������� �������������

	double *curDist;	// ��������������� ������ ��� ���������

	unsigned int *nu;	// ������� ������ ��� �������� ���������� �������������
						// � ����������� �� ������� ������

	double get_buffer_likelihood(Buffer &b);
	/* ����������� ������� ���������� ��� ����������� ������������� */
	double fix_computational_error(double *arr, int size);

	/* ��� �������� nu, ��� ��. ������� T � ����� s*/
	double vectorIteration(double *target, double *deriv, int size,
		double eps, double &step, double lh0);
	/* likelihood ����������� �� ������������������ */
	double vectorIteration(double *target, double *deriv, int size,
		double eps, double &step, double lh0, istream *is, int s_pos, int e_pos);

	/*derivatives !!! �� ��� ������� ��� ������, ������ Q 
	(� ���� ����� - ���� �� ����)*/
	void lambda_derivatives(double *dLambda);	// �� nu
	void lambda_derivatives(double *dLambda, istream *is, int s_pos, int e_pos);

	/* ����������� ������ ������� Q � ������� l*/
	void Q_derivatives(double *dQl, int l);		// �� nu
	void Q_derivatives(double *dQl, int l, istream *is, int s_pos, int e_pos);

	/* � ������� nu */
	double iterativeEstimation(ostream &os, double eps = 0.001, int max_iter = 100,
		double start_step = 0.1, bool report = 0);

public:
	MTDModel();
	MTDModel(char *fname);					// ��� ���������
	MTDModel(int s, int L, char *alpha);	// ��� ������ ����������, � ���������� �� ���������
	~MTDModel();

	void estimateNu(istream *is);
	void estimateNu(istream *is, int s_pos, int e_pos);

	void estimateInitialParameters(istream *is);
	void estimateInitialParameters(istream *is, int s_pos, int e_pos);

	/* � ������� nu */
	double iterativeEstimation_Nu(istream *is, int s_pos, int e_pos, ostream &os,
		double eps = 0.001, int max_iter = 100,	double start_step = 0.1, bool report = 0);
	/* ��� ������ nu */
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
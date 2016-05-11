#include "MCModel.h"

#ifndef JLMODEL
#define JLMODEL
/*
������ �������� �����:
������� ������
����� ���������
�������
������������� �� ����� ������
������������� ������ ����� ������
����� ��
*/
class JacobsLewisModel : public MarkovChainModel
{
	double ro;			// �������� ro
	double *lambda;		// ������ ������������ lambda
	double *P;			// ������ ������������ pi

	// ������������ ���� (��� ������� ���������� ������)
	unsigned int *nu;			// ������ ��� �������� ������� ������

	void estimateNu(istream *is); 
	double *estimateQ(double *Q); // ����� ������ Nu

	double estimateInitialRo(double *Q); // ����� ������ P
	double *estimateInitialLambda(double *Q); // ����� ����� P

	double get_buffer_likelihood(Buffer &b);// ��������������� ��� ��������� 6
	void P_derivatives(double* d_P);	// ����� ������ nu, P, lambda, ro
	void Lambda_derivatives(double* d_lambda);
	double Ro_derivative();
	void all_derivatives(double *d_P, double *d_lambda, double &d_ro); 

	void P_derivatives(double* d_P, istream *is, int s_pos, int e_pos);
	void Lambda_derivatives(double* d_lambda, istream *is, int s_pos, int e_pos);
	double Ro_derivative(istream *is, int s_pos, int e_pos);
	void all_derivatives(double *d_P, double *d_lambda, double &d_ro, 
		istream *is, int s_pos, int e_pos); // ���� ����� ��������� �� �� ���� ������

	/*one iteration for vectorized parameter, returns new likelihood*/
	double vectorIteration(double *target, double *deriv, int size, 
		double eps, double &step, double lh0);
	/*one iteration for single parameter, returns new likelihood*/
	double singleIteration(double & target, double deriv,
		double eps, double &step, double lh0);
	
	/*one iteration for vectorized parameter, returns new likelihood*/
	double vectorIteration(double *target, double *deriv, int size,
		double eps, double &step, double lh0, istream *is, int s_pos, int e_pos);
	/*one iteration for single parameter, returns new likelihood*/
	double singleIteration(double & target, double deriv,
		double eps, double &step, double lh0, istream *is, int s_pos, int e_pos);

	__int64 next1(Buffer &b, int i); // s+1-� ���������� ������ ���� ����� i, ��������� != i
	__int64 next2(Buffer &b, int i); // ��� ���������� ����� i

public:
	JacobsLewisModel();
	JacobsLewisModel(int s, int L, char* alpha);	// ������ ������, ��� ������ ����������;
													// lambda = (1/s), P = (1/L), ro = 0.5
	JacobsLewisModel(char *fname);					// ���������������� ������ ����������� �� �����
	~JacobsLewisModel();

	void estimateP(istream *is);
	void estimateInitialParameters(istream *is);	// !����� ������ nu � P

	/* ��� ������������ ����������*/
	void iterativeEstimation(istream *is, ostream &os, double eps = 0.001);
	/* ��� ������������� ���������� */
	double iterativeEstimation(istream *is, int s_pos, int e_pos,
		double eps = 0.001, int max_iter = 40);
	/* �������� ������, �� ������������ �������� � 3 ���� ������, � �������*/
	void iterativeEstimation(istream *is, int s_pos, int e_pos,
		ostream &os, double eps = 0.001);

	double likelihood();	// �� ������� ������
	double likelihood(istream *is);
	double likelihood(istream *is, int s_pos, int e_pos);
	__int64 numberOfParams();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
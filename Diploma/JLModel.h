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

	double get_buffer_likelihood(Buffer b);// ��������������� ��� ��������� 6
	void P_derivatives(double* d_P);	// ����� ������ nu, P, lambda, ro
	void Lambda_derivatives(double* d_lambda);
	double Ro_derivative();
	void all_derivatives(double *d_P, double *d_lambda, double &d_ro); 

	void P_derivatives(double* d_P, istream *is, int s_pos, int e_pos);
	void Lambda_derivatives(double* d_lambda, istream *is, int s_pos, int e_pos);
	double Ro_derivative(istream *is, int s_pos, int e_pos);
	void all_derivatives(double *d_P, double *d_lambda, double &d_ro, 
		istream *is, int s_pos, int e_pos); // ���� ����� ��������� �� �� ���� ������

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

	double likelihood(istream *is);
	__int64 numberOfParams();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
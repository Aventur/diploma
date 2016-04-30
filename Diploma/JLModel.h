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

	double likelihood();
	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
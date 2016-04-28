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

	unsigned int *nu;			// ������ ��� �������� ������� ������

public:
	JacobsLewisModel();
	JacobsLewisModel(int s, int L, char* alpha);	// ������ ������, ��� ������ ����������;
													// lambda = (1/s), P = (1/L), ro = 0.5
	JacobsLewisModel(char *fname);					// ���������������� ������ ����������� �� �����
	~JacobsLewisModel();

	void estimateNu(istream &is);
	void estimateP(istream &is);

	double likelihood();
	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
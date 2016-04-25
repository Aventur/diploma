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

	int *nu;			// ������ ��� �������� ������� ������

public:
	JacobsLewisModel();
	JacobsLewisModel(int s, int L, char* alpha);	// ������ ������, ��� ������ ����������;
													// lambda = (1/s), P = (1/L), ro = 0.5
	JacobsLewisModel(char *fname);					// ���������������� ������ ����������� �� �����
	~JacobsLewisModel();



	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
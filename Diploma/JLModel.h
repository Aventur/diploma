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

public:
	JacobsLewisModel();
	JacobsLewisModel(char *fname);
	~JacobsLewisModel();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
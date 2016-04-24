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

	double *curDist;	// ��������������� ������

public:
	MTDModel();
	MTDModel(char *fname);
	~MTDModel();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
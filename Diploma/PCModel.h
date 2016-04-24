#include "MCModel.h"

#ifndef PCMODEL
#define PCMODEL
/*
������ �������� �����: (��������������� �������� full)
������� ������ (0, 1, 2)
����� ��������� (all)
������� (all)
������������� �� ����� ������ (0)
����� ������ r (0, 1, 2)
������ ������ (r ����� �� 0..s-1) (0, 1)
������� Q �� ������� �� L ����� (0)
*/
class PartialConnectionsModel : public MarkovChainModel
{
	double *P;			// ���������� ������������ �������������
	double **Q;			// ������� ��������� � ��������� �����
	int r;				// ����� ������
	int *m;				// ������ ������
	int c;				// ����� ����� � ������� Q

	int number();		// ���������� ����� ������� ��������������������� � ������ � ������ �������
	bool nextM(int level);

public:
	PartialConnectionsModel();
	PartialConnectionsModel(char *fname, short full = 0);
	~PartialConnectionsModel();

	void setParam(int p1, int p2);
	void estimateModel(istream &is);	// �������� ������ s, r, L, alpha (full = 2)
	double estimateModelQ(istream &is);	// �������� ��, ����� Q (full = 1)
	int nextInitialState();
	int nextState();
	void printModel(ostream &os, short full);
};

#endif
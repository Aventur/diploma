#include "MCModel.h"

#ifndef COMODEL
#define COMODEL
/*
������ �������� �����:
������� ������ (0, 1, 2)
����� ��������� (all)
������� (all)
������������� �� ����� ������ (0)
����� ��� (���� �����) (0, 1, 2)
����� M ��������� ������ Q (���� �����) (0, 1)
������ ������, ��������������� ��� (K ����� mk) (0, 1)
������ ���������, �� ������� �������������� �������, ��� ������� ��� (K ����� bk) (0, 1)
M ������ Q �� ������� (0)
*/
class ConditionalOrderModel : public MarkovChainModel
{
	double *P;			// ������������ �������������
	int B;				// ����� ���
	int M;				// ����� ��������� ������ Q
	int K;				// ����� ��������� ��������� ���
	int *mk;			// ������ ������, ��������������� ��� (������� � 0)
	int *bk;			// ������ ���������, �� ������� �������������� �������, �����. ��� (������� � 0)
	double **Q;			// M ������ � �����, �� �������

	int number();		// ����� �������� ��� (����, ��� � ������)

public:
	ConditionalOrderModel();
	ConditionalOrderModel(char *fname, short full = 0);
	~ConditionalOrderModel();

	void setParam(int p1, int p2);
	void estimateModel(istream &is);
	void estimateModelB(istream &is);
	void estimateModelQ(istream &is);	// ��� �������� alphabet, L, s, B, M, mk, bk
	int nextInitialState();
	int nextState();
	void printModel(ostream &os, short full = 0);
};

#endif
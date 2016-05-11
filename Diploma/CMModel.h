#include "MCModel.h"

#ifndef CLASSICMARKOVMODEL
#define CLASSICMARKOVMODEL
/*
������ �������� �����:
������� ������
����� ���������
�������
������� P �� �������
*/
class ClassicMarkovModel : public MarkovChainModel
{
public:
	double *Q;					// ������� ������������ ���������, ���������� � ������
	unsigned int *nu;			// ������ ��� �������� ������� ������

	void estimateNu(istream *is);
	void estimateNu(istream *is, int s_pos, int e_pos);

	double get_buffer_likelihood(Buffer &b);// ���������������

//public:
	ClassicMarkovModel();
	ClassicMarkovModel(int s, int L, char* alpha);	// ������ ������, ��� ������ ����������;
	ClassicMarkovModel(char *fname);	// ���������������� ������ ����������� �� �����
	~ClassicMarkovModel();

	double estimateQ(istream *is);
	/*returns likelihood of the estimated model*/
	double estimateQ(istream *is, int s_pos, int e_pos);

	double likelihood(istream *is);
	double likelihood(istream *is, int s_pos, int e_pos);
	__int64 numberOfParams();

	int nextInitialState();
	int nextState();
	void printModel(ostream &os);
};

#endif
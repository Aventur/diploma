#include "bsv.h"
#include <fstream>
#include <sstream>
#include <time.h>

#ifndef MARKOVCHAINMODEL
#define MARKOVCHAINMODEL

class MarkovChainModel
{
protected:
	BSV *bsv;			// ��������� ���
	int s;				// ������� ������
	int L;				// ���������� ���������
	char *alphabet;		// ��������� ��������, ��������������� ����������
	int *buffer;		// ��������� ��������������� s ���������

	void shiftBuffer(int state);

public:
	MarkovChainModel();
	~MarkovChainModel();

	int indexof(char t);	// ����� ������� � ��������
	virtual void setParam(int p1, int p2);
	virtual void estimateModel(istream &is);	// ����� ���������� L, alphabet
	virtual int nextInitialState() = 0;
	virtual int nextState() = 0;
	virtual void printModel(ostream &os, short full = 0);
	void generateSequence(int n, ostream &os);
	static void printArray(double *arr, int n, ostream &os);
	static void printArray(int *arr, int n, ostream &os);
	int nextDiscreteRand(double *dist, int n);
};

#endif
#include "bsv.h"
#include "mcutils.h"
#include <fstream>
#include <sstream>
#include <time.h>
#include <Windows.h>

#ifndef MARKOVCHAINMODEL
#define MARKOVCHAINMODEL

struct Buffer {
	int *str;
	int L;
	int len;
	int beg;
	int end;

	Buffer(int L, int len);
	~Buffer();

	__int64 getNumber();
	void shiftBuffer(int i);
};

class MarkovChainModel
{
protected:
	BSV *bsv;			// ��������� ���
	int s;				// ������� ������
	int L;				// ���������� ���������
	char *alphabet;		// ��������� ��������, ��������������� ����������
	int *buffer;		// ��������� ��������������� s ���������

	void shiftBuffer(int state);
	int indexof(char t);// ������� ��� ����������� ������ ������ P, �����.
	int number(char* t);// �������� ������������������ ���������

public:
	MarkovChainModel();
	~MarkovChainModel();

	virtual double likelihood() = 0;
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
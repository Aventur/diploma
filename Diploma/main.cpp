#include "MCModel.h"
#include "COModel.h"
#include <fstream>
#include <conio.h>

int main()
{
	/*
	MarkovChainModel *co;
	co = new ConditionalOrderModel("COinput.txt");

	ofstream os("COoutput5.txt");
	co->generateSequence(100000, os);
	co->printModel(cout);
	os.close();

	ifstream is("COoutput5.txt");
	co->setParam(12, 1);
	co->estimateModel(is);
	is.close();

	os.open("COo5_res.txt");
	co->printModel(os);
	os.close();

	if (co != NULL) delete co; */

	Buffer b = Buffer(4, 10);
	b.shiftBuffer(2);
	b.shiftBuffer(1);

	for (int i = 0; i < 10; i++)
		cout << b.str[i];

	cout << '\n' << b.getNumber() << endl;
	_getch();
}


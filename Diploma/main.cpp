#include "MCModel.h"
#include "COModel.h"
#include "JLModel.h"
#include "mcutils.h"
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

	/*
	JacobsLewisModel *m = new JacobsLewisModel("JLinput.txt");

	ofstream f;
	f.open("genJL.txt");

	m->generateSequence(10000000, f);

	f.close();

	delete m;
	*/

	
	fstream is;
	is.open("genJL.txt", ios_base::in);

	JacobsLewisModel *jl = new JacobsLewisModel(5, 4, "ACGT");
	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	jl->estimateInitialParameters(sis);

	jl->printModel(cout);

	_getch();

	delete jl;
	
	_getch();
}


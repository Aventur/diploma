#include "MCModel.h"
#include "COModel.h"
#include "JLModel.h"
#include "MTDModel.h"
#include "mcutils.h"
#include <fstream>
#include <conio.h>
#include <time.h>


double *JL_KOP_exp(char *filename, ostream &os, int s, int bulk_size, int step, int n_steps=100)
{
	fstream is;
	is.open(filename, ios_base::in);

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	istringstream *is1, *is2, *is0;
	double l1, l2, l0;
	JacobsLewisModel *jl = new JacobsLewisModel(s, 4, "ACGT");

	double *KOP = new double[n_steps];
	for (int i = 0; i < n_steps * step; i += step)
	{
		cout << ((!(i % (step * 100))) ? (i / step) : 0 );	//logging
		os << "beg_pos " << i << ", bulk_size = " << bulk_size << endl;
		is1 = new istringstream(sis->str().substr(i, bulk_size));
		is2 = new istringstream(sis->str().substr(i + bulk_size, bulk_size));
		is0 = new istringstream(sis->str().substr(i, 2 * bulk_size));

		jl->estimateInitialParameters(is1);
		jl->printModel(os);
		l1 = jl->likelihood(is1);
		jl->estimateInitialParameters(is2);
		jl->printModel(os);
		l2 = jl->likelihood(is2);
		jl->estimateInitialParameters(is0);
		jl->printModel(os);
		l0 = jl->likelihood(is0);

		if (is1 != NULL) delete is1;
		if (is2 != NULL) delete is2;
		if (is0 != NULL) delete is0;

		KOP[i / step] = 2. * (l1 + l2 - l0);
		os << l1 << " " << l2 << " " << l0 << endl;
		os << "KOP value" << " " << KOP[i/step] << endl;
		os << "*********************************************" << endl;
	}
	cout << "\n\n\n\n";

	for (int i = 0; i < n_steps; i++)
		os << KOP[i] << " ";
	os << "\nr = " << s + 3 << endl;
	return KOP;
}

double *JL_KOP_exp_stream(char *filename, int s_pos, ostream &os, int s, 
	int bulk_size, int step, int n_steps = 100, double eps = 0.001, int max_iter = 40)
{
	fstream is;
	is.open(filename, ios_base::in);

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	double l1, l2, l0;
	JacobsLewisModel *jl1 = new JacobsLewisModel(s, 4, "ACGT");
	JacobsLewisModel *jl2 = new JacobsLewisModel(s, 4, "ACGT");
	JacobsLewisModel *jl0 = new JacobsLewisModel(s, 4, "ACGT");

	double *KOP = new double[n_steps];
	for (int i = s_pos; i < s_pos + n_steps * step; i += step)
	{
		cout << ((!((i - s_pos) % (step * 100))) ? ((i  - s_pos)/ step) : 0);//logging
		os << "beg_pos " << i << ", bulk_size = " << bulk_size << endl;

		l1 = jl1->iterativeEstimation(sis, i, i + bulk_size, eps, max_iter);
		//jl->printModel(os);

		l2 = jl2->iterativeEstimation(sis, i + bulk_size, i + 2 * bulk_size, eps, max_iter);
		//jl->printModel(os);

		l0 = jl0->iterativeEstimation(sis, i, i + 2 * bulk_size, eps, max_iter);
		//jl->printModel(os);

		KOP[(i - s_pos) / step] = 2. * (l1 + l2 - l0);
		os << l1 << " " << l2 << " " << l0 << endl;
		os << "KOP value" << " " << KOP[(i - s_pos) / step] << endl;
		os << "*********************************************" << endl;
	}
	cout << "\n\n\n\n";

	for (int i = 0; i < n_steps; i++)
		os << KOP[i] << " ";
	os << "\nr = " << s + 3 << endl;
	if (jl1 != NULL) delete jl1;
	if (jl2 != NULL) delete jl2;
	if (jl0 != NULL) delete jl0;
	return KOP;
}

void JL_analyze(char *filename, ostream &os, int MAX_S = 12)
{
	fstream is;
	is.open(filename, ios_base::in);

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	JacobsLewisModel **jls = new JacobsLewisModel*[MAX_S];
	for (int i = 0; i < MAX_S; i++)
		jls[i] = new JacobsLewisModel(i + 1, 4, "ACGT");

	double mAIC, mBIC;
	int argminAIC, argminBIC;
	mAIC = DBL_MAX; mBIC = DBL_MAX;
	argminAIC = -1; argminBIC = -1;

	double aic, bic;
	for (int i = 0; i < MAX_S; i++)
	{
		jls[i]->estimateInitialParameters(sis);
		jls[i]->printModel(os);
		aic = jls[i]->AIC(sis);
		bic = jls[i]->BIC(sis);

		os << "AIC: " << aic << endl;
		os << "BIC: " << bic << endl;
		os << "*****************************************" << endl;

		if (aic < mAIC) { mAIC = aic; argminAIC = i + 1; }
		if (bic < mBIC) { mBIC = bic; argminBIC = i + 1; }
	}

	os << "argminAIC: " << argminAIC << endl;
	os << "argminBIC: " << argminBIC << endl;

	for (int i = 0; i < MAX_S; i++)
		if (jls[i] != 0) delete jls[i];
	if (jls != NULL) delete[]jls;
}






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
	f.open("genJL2.txt");

	m->generateSequence(30000, f);
	delete m;
	m = new JacobsLewisModel("JLinput2.txt");
	m->generateSequence(30000, f);

	f.close();

	delete m;
	*/
	
	/*
	JacobsLewisModel *jl = new JacobsLewisModel(5, 4, "ACGT");
	jl->printModel(cout);
	jl->estimateInitialParameters(sis);

	jl->printModel(cout);

	_getch();

	delete jl;*/
	//JL_analyze("genJL.txt", cout);
	/*
	ofstream os;
	double *k;
	
	os.open("KOP_genJL.txt");
	k = JL_KOP_exp("genJL.txt", os, 5, 10000, 3, 15000);
	if (k != NULL) delete[]k;
	os.close();
	
	os.open("KOP_genJL2.txt");	//разладка на 30 005 эл-те
	k = JL_KOP_exp("genJL2.txt", os, 5, 10000, 3, 15000);
	if (k != NULL) delete[]k;
	os.close();

	os.open("KOP_clx5.txt");
	k = JL_KOP_exp("callithrix.txt", os, 5, 10000, 3, 15000);
	if (k != NULL) delete[]k;
	os.close();

	os.open("KOP_clx10.txt");
	k = JL_KOP_exp("callithrix.txt", os, 10, 10000, 3, 15000);
	if (k != NULL) delete[]k;
	os.close();*/

	/*
	int L = 4;
	int s = 5;

	ofstream os;
	double *k;
	
	os.open("KOP_clx5_110-140k_s30_w500.txt");
	k = JL_KOP_exp_stream("callithrix.txt", 110000, os, s, 500, 30, 1000, 0.0005, 80);
	if (k != NULL) delete[]k;
	os.close();

	os.open("KOP_clx5_435-465k_s30_w500.txt");
	k = JL_KOP_exp_stream("callithrix.txt", 435000, os, s, 500, 30, 1000, 0.0005, 80);
	if (k != NULL) delete[]k;
	os.close();

	os.open("KOP_clx5_1120-1150k_s30_w500.txt");
	k = JL_KOP_exp_stream("callithrix.txt", 1120000, os, s, 500, 30, 1000, 0.0005, 80);
	if (k != NULL) delete[]k;
	os.close();

	cout << "STOP!!!!!" << endl;*/

	/* reading from file, end position
	ifstream is("MTDgen.txt");
	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();


	sis->seekg(0, ios::end);
	int T = sis->tellg();
	cout << T << endl;
	_getch();
	sis->seekg(0, ios::beg);

	for (int i = 0; i <= T; i++)
		cout << sis->get() << endl;/**/


	MTDModel mtd_gen("MTD_input.txt");

	ofstream os("MTDgen.txt");
	mtd_gen.generateSequence(30000, os);
	os.close();

	MTDModel mtd_est(5, 4, "ACGT");
	ifstream is("MTDgen.txt");

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	cout << mtd_est.likelihood(sis, 0, 30000) << endl;
	mtd_est.estimateInitialParameters(sis, 0, 30000);
	mtd_est.printModel(cout);
	cout << mtd_est.likelihood(sis, 0, 30000) << endl;

	if (sis != NULL) delete sis;/**/
	_getch();
}


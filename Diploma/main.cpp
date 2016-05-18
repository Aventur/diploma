#include "MCModel.h"
#include "COModel.h"
#include "PCModel.h"
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

void JL_analyze(char *filename, int s_pos, int e_pos, ostream &os, int MAX_S = 12)
{
	fstream is;
	is.open(filename, ios_base::in);

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	JacobsLewisModel **jls = new JacobsLewisModel*[MAX_S];
	for (int i = 0; i < MAX_S; i++)
		jls[i] = new JacobsLewisModel(i + 1, 4, "ACGT");

	double mAIC, mBIC;
	int argminAIC, argminBIC, s_time;
	mAIC = DBL_MAX; mBIC = DBL_MAX;
	argminAIC = -1; argminBIC = -1;

	os << "JacobsLewisModel analysis. File " << filename << "; s_pos = " <<
		s_pos << ", e_pos = " << e_pos << ", len = " << e_pos - s_pos << ".\n";
	double aic, bic;
	for (int i = 0; i < MAX_S; i++)
	{
		cout << i << endl;
		s_time = time(NULL);
		if (i <= 6)
		{
			//jls[i]->estimateInitialParameters(sis);
			jls[i]->iterativeEstimation(sis, s_pos, e_pos, 0.00005, (i+2) * 15);
		}
		else
			jls[i]->iterativeEstimation(sis, s_pos, e_pos, 0.00005, (i + 2) * 15);


		jls[i]->printModel(os);
		aic = jls[i]->AIC(sis, s_pos, e_pos);
		bic = jls[i]->BIC(sis, s_pos, e_pos);

		os << "AIC: " << aic << endl;
		os << "BIC: " << bic << endl;
		os << "Time to analyze: " << time(NULL) - s_time << " s.\n";
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

void JL_analyze_with_random(char *filename, int s_pos, int e_pos, ostream &os, int MAX_S = 12, int rnd_inits=10)
{
	fstream is;
	is.open(filename, ios_base::in);

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	JacobsLewisModel **jls = new JacobsLewisModel*[MAX_S];
	for (int i = 0; i < MAX_S; i++)
		jls[i] = new JacobsLewisModel(i + 1, 4, "ACGT");

	double mAIC, mBIC;
	int argminAIC, argminBIC, s_time;
	mAIC = DBL_MAX; mBIC = DBL_MAX;
	argminAIC = -1; argminBIC = -1;

	os << "JacobsLewisModel analysis. File " << filename << "; s_pos = " <<
		s_pos << ", e_pos = " << e_pos << ", len = " << e_pos - s_pos << ".\n";
	double aic, bic;
	for (int i = 0; i < MAX_S; i++)
	{
		cout << i << endl;
		s_time = time(NULL);
		JacobsLewisModel *jl_exp = new JacobsLewisModel(i + 1, 4, "ACGT");
		JacobsLewisModel *jl_swap;
		double lh, mlh = -DBL_MAX;

		mlh = jls[i]->iterativeEstimation(sis, s_pos, e_pos, 0.0005, (i + 2) * 20);
		for (int j = 0; j < rnd_inits; j++)
		{
			cout << j;
			jl_exp->setRandomInitialParameters();
			lh = jl_exp->iterativeEstimation(sis, s_pos, e_pos, 0.0005, (i + 2) * 20);
			if (lh > mlh) { jl_swap = jl_exp; jl_exp = jls[i]; jls[i] = jl_swap; mlh = lh; }
		}
		if (jl_exp != NULL) delete jl_exp;


		jls[i]->printModel(os);
		aic = jls[i]->AIC(sis, s_pos, e_pos);
		bic = jls[i]->BIC(sis, s_pos, e_pos);

		os << "AIC: " << aic << endl;
		os << "BIC: " << bic << endl;
		os << "Time to analyze: " << time(NULL) - s_time << " s.\n";
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

void MTD_analyze(char *filename, int s_pos, int e_pos, ostream &os, int MAX_S = 20)
{
	fstream is;
	is.open(filename, ios_base::in);

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	MTDModel **mtds = new MTDModel*[MAX_S];
	for (int i = 0; i < MAX_S; i++)
		mtds[i] = new MTDModel(i + 1, 4, "ACGT");

	double mAIC, mBIC;
	int argminAIC, argminBIC, s_time;
	mAIC = DBL_MAX; mBIC = DBL_MAX;
	argminAIC = -1; argminBIC = -1;

	double aic, bic;
	os << "MTDModel analysis. Source file " << filename << "; s_pos = " <<
		s_pos << ", e_pos = " << e_pos << ", len = " << e_pos - s_pos << ".\n";
	for (int i = 0; i < MAX_S; i++)
	{
		cout << i << endl;
		s_time = time(NULL);
		mtds[i]->estimateInitialParameters(sis, s_pos, e_pos);
		
		if (pow(4, i + 2) > e_pos - s_pos)
			mtds[i]->iterativeEstimation(sis, s_pos, e_pos, os, 0.00005, (i+2)*10, 0.02, false);
		else
			mtds[i]->iterativeEstimation_Nu(sis, s_pos, e_pos, os, 0.00005, (i+2)*10, 0.02, false);
		
		mtds[i]->printModel(os);
		aic = mtds[i]->AIC(sis, s_pos, e_pos);
		bic = mtds[i]->BIC(sis, s_pos, e_pos);

		os << "AIC: " << aic << endl;
		os << "BIC: " << bic << endl;
		os << "Time to analyze: " << time(NULL) - s_time << " s.\n";
		os << "*****************************************" << endl;

		if (aic < mAIC) { mAIC = aic; argminAIC = i + 1; }
		if (bic < mBIC) { mBIC = bic; argminBIC = i + 1; }
	}

	os << "argminAIC: " << argminAIC << endl;
	os << "argminBIC: " << argminBIC << endl;

	for (int i = 0; i < MAX_S; i++)
		if (mtds[i] != 0) delete mtds[i];
	if (mtds != NULL) delete[]mtds;
}

void PC_analyze(char *filename, ostream &os, int MAX_S = 40, int MAX_R = 3)
{
	fstream is;
	is.open(filename, ios_base::in);

	//istringstream *sis = get_string_stream(is, "ACGT");
	//is.close();

	PartialConnectionsModel pc("PC_input.txt", 2);
	int s_time = 0;

	//for (int i = 0; i < MAX_R; i++)
	//{
		//cout << i << endl;
		s_time = time(NULL);
		//pc.setParam(MAX_S, i + 1);
		pc.estimateModel(is);

		pc.printModel(os, 1);
		os << "Time to analyze: " << time(NULL) - s_time << " s.\n";
		os << "*****************************************" << endl;
	//}
		//if (sis != NULL) delete sis;
}

void exp_JL_MTD_3species()
{
	/* JL analyze /*/
	ofstream os("JL_analyze_clx.txt");
	JL_analyze("clx.txt", 0, 10000, os, 25);
	os.close();	/**/
				/* JL analyze /*/
	os.open("JL_analyze_cns.txt");
	JL_analyze("cns.txt", 0, 10000, os, 25);
	os.close();	/**/
				/* JL analyze /*/
	os.open("JL_analyze_hsp.txt");
	JL_analyze("hsp.txt", 0, 10000, os, 25);
	os.close();	/**/

	/* MTD analyze /*/
	os.open("MTD_analyze_clx.txt");
	MTD_analyze("clx.txt", 0, 10000, os, 20);
	os.close();	/**/
				/* MTD analyze /*/
	os.open("MTD_analyze_cns.txt");
	MTD_analyze("cns.txt", 0, 10000, os, 20);
	os.close();	/**/
				/* MTD analyze /*/
	os.open("MTD_analyze_hsp.txt");
	MTD_analyze("hsp.txt", 0, 10000, os, 20);
	os.close();	/**/
}

void exp_PC_3species()
{
	/* PC analyze /*/
	ofstream os("PC_analyze_clx.txt");
	PC_analyze("clx.txt", os, 20, 3);
	os.close();	/**/
				/* PC analyze /*/
	os.open("PC_analyze_cns.txt");
	PC_analyze("cns.txt", os, 20, 3);
	os.close();	/**/
				/* PC analyze /*/
	os.open("PC_analyze_hsp.txt");
	PC_analyze("hsp.txt", os, 20, 3);
	os.close();	/**/
}

void exp_reest_callithrix()
{
	/* JL analyze /*/
	ofstream os("JL_analyze_callithrix.txt");
	JL_analyze("callithrix.txt", 0, 200000, os, 20);
	os.close();	/**/
				/* MTD analyze /*/
	os.open("MTD_analyze_callithrix.txt");
	MTD_analyze("callithrix.txt", 0, 200000, os, 20);
	os.close();	/**/
}

void exp_generated_jl_mtd()
{
	/* MTD analyze /*/
	ofstream os("MTD_analyze_gen.txt");
	MTD_analyze("MTDgen.txt", 0, 100000, os, 30);
	os.close();	/**/

	/* JL analyze /*/
	os.open("JL_analyze_gen.txt");
	JL_analyze("JLgen.txt", 0, 100000, os, 30);
	os.close();	/**/
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
	/* JL KOP experiments
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

	/* JL KOP experiments
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
	/**/

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

	/* MTD estimation testing /*
	MTDModel mtd_est(6, 4, "ACGT");
	ifstream is("MTDgen.txt");

	istringstream *sis = get_string_stream(is, "ACGT");
	is.close();

	cout << mtd_est.likelihood(sis, 0, 30000) << endl;

	mtd_est.estimateInitialParameters(sis);
	mtd_est.printModel(cout);
	cout << mtd_est.likelihood(sis, 0, 30000) << endl;

	ofstream os("MTD_est_report.txt");
	int s_time = time(NULL);
	mtd_est.iterativeEstimation(sis, 0, 30000, os, 0.0001, 100, 0.02, true);
	mtd_est.printModel(cout);
	cout << mtd_est.likelihood(sis, 0, 30000) << endl;
	cout << time(NULL) - s_time << endl;
	os.close();

	if (sis != NULL) delete sis;/**/

	/* MTD generation /*

	int T = 100000;
	MTDModel mtd_gen("MTD_input.txt");

	ofstream os("MTDgen.txt");
	mtd_gen.generateSequence(T, os);
	os.close();/**/

	/* JL generation /*

	int T = 100000;
	JacobsLewisModel jl_gen("JL_input.txt");

	ofstream os("JLgen.txt");
	jl_gen.generateSequence(T, os);
	os.close();/**/

	/* MTD analyze /*
	ofstream os("MTD_analyze_gen.txt");
	MTD_analyze("MTDgen.txt", 0, 100000, os, 40);
	os.close();	/**/

	/* JL analyze /*
	ofstream os("JL_analyze_gen.txt");
	JL_analyze("JLgen.txt", 0, 100000, os, 40);
	os.close();	/**/
	
	//exp_JL_MTD_3species();
	//cout << "----------___-------_____--------------------_--__-__------" << endl;
	//exp_generated_jl_mtd();
	//exp_PC_3species();

	ofstream os("jl_exp_randominit_clx.txt");
	//JL_analyze_with_random("clx.txt", 0, 100000, os, 20, 10);
	os.close();

	cout << "STOP!!!!!" << endl;
	_getch();
}


	/* main for generating all types of sequences (testing)
	MarkovChainModel *mtd;
	MarkovChainModel *jl;
	MarkovChainModel *pc;
	MarkovChainModel *co;

	ifstream is("Streptococcus.txt");

	mtd = new MTDModel("MTDinput.txt");
	mtd->printModel(cout);
	mtd->generateSequence(1000, cout);
	cout << endl;
	
	jl = new JacobsLewisModel("JLinput.txt");
	jl->printModel(cout);
	jl->generateSequence(1000, cout);
	cout << endl;

	pc = new PartialConnectionsModel("PCinput.txt");
	pc->printModel(cout);
	pc->generateSequence(1000, cout);
	cout << endl;
	pc->estimateModel(is);
	pc->printModel(cout);


	co = new ConditionalOrderModel("COinput.txt");
	co->printModel(cout);
	co->generateSequence(1000, cout);
	cout << endl;

	getch();
	is.close();
	if (mtd != NULL)
		delete mtd;
	if (jl != NULL)
		delete jl;
	if (pc != NULL)
		delete pc;
	if (co != NULL)
		delete co;*/


/*	тестирование генерации и оценки без пор€дка дл€ PCModel
	MarkovChainModel *pc;
	pc = new PartialConnectionsModel("PCinput.txt");

	ofstream os("PCoutput5.txt");
	pc->generateSequence(100000, os);
	pc->printModel(cout);
	os.close();

	ifstream is("PCoutput5.txt");
	pc->estimateModel(is);
	is.close();
	pc->printModel(cout);

	if (pc != NULL) delete pc;
	getch();

	*/

//Homo sapiens 1 chrm, NW_004929289.1


/*		CO_experiments.txt generator ;)
	MarkovChainModel *co;
	co = new ConditionalOrderModel("COinput1.txt", 2);
	
	ifstream is("callithrix.txt");
	ofstream os("CO_experiments.txt");
	os << "Callithrix.txt, NW_003182879.1\n";
	int s[] = {12, 12, 12, 12, 32, 32, 32, 32, 64, 64, 64, 64};
	int b[] = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};

	for (int i = 0; i < 12; i++)
	{
		co->setParam(s[i], b[i]);
		co->estimateModel(is);
		is.clear();
		is.seekg(0, ios::beg);
		co->printModel(os, -1);
		cout << "+";
	}
	is.close();
	is.open("Streptococcus.txt");
	os << "\n\nStreptococcus.txt, whole sequence NC_021900.1\n";
	cout << endl;

	for (int i = 0; i < 12; i++)
	{
		co->setParam(s[i], b[i]);
		co->estimateModel(is);
		is.clear();
		is.seekg(0, ios::beg);
		co->printModel(os, -1);
		cout << "+";
	}
	is.close();
	is.open("COoutput5.txt");
	os << "\n\nCOoutput5.txt, generated sequence, s=10, B=1, bk=7 7 5 5\n";
	cout << endl;

	for (int i = 0; i < 12; i++)
	{
		co->setParam(s[i], b[i]);
		co->estimateModel(is);
		is.clear();
		is.seekg(0, ios::beg);
		co->printModel(os, -1);
		cout << "+";
	}
	cout << endl;

	is.close();
	os.close();
	if (co != NULL) delete co;
	getch();
*/
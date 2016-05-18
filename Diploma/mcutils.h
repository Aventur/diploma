#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

#ifndef MARKOVCHAINSUTILS
#define MARKOVCHAINSUTILS

template <class T> int cmp(const void *a, const void *b)
{
	if (*(T*)a < *(T*)b) return -1;
	if (*(T*)a > *(T*)b) return 1;
	if (*(T*)a == *(T*)b) return 0;
}

istringstream *get_string_stream(istream &is, const char* alphabet);
istream *refresh_string_stream(istream *is);
void project_to_simplex(double *vector, int size);	// только нормированный в L1 вектор!





#endif
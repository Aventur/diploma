#include <fstream>
#include <sstream>
using namespace std;

#ifndef MARKOVCHAINSUTILS
#define MARKOVCHAINSUTILS

istringstream *get_string_stream(istream &is, const char* alphabet);
istream *refresh_string_stream(istream *is);






#endif
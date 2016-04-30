#include "mcutils.h"

istringstream * get_string_stream(istream &is, const char * alphabet)
{
	is.seekg(0, ios::end);						// перенос последовательности
	int g = is.tellg();							// из файла в оперативную
	is.seekg(0, ios::beg);						// память
	char *file = new char[g];
	is.read(file, g);

	int new_g = 0;								// считаем количество допустимых символов
	for (int i = 0; i < g; i++)
		if (strchr(alphabet, file[i]) != NULL)
			new_g++;
	char *new_file = new char[new_g];			// новая строка, для допустимых символов

	new_g = 0;
	for (int i = 0; i < g; i++)					// перенос допустимых символов
		if (strchr(alphabet, file[i]) != NULL)
			new_file[new_g++] = file[i];

	istringstream *sis = new istringstream(new_file);
	delete[]file;
	delete[]new_file;

	return sis;
}

istream * refresh_string_stream(istream *is)
{
	is->clear(0);
	is->seekg(0, ios::beg);
	return is;
}

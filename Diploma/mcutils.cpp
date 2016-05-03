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

void project_to_simplex(double *vector, int size)
{
	int i, j;
	int *order = new int[size];
	for (int i = 0; i < size; i++)
		order[i] = i;

	int t1;
	double sum;					// сортируем вставками, запоминая порядок
	for (int i = 1; i < size; i++)
	{
		for (j = i - 1; vector[i] > vector[j] && j >= 0; j--);
		if (j == -1) continue;
		else
		{
			t1 = order[i]; order[i] = order[j]; order[j] = t1;
			sum = vector[i]; vector[i] = vector[j]; vector[j] = sum;
		}
	}

	i = 0; sum = 0;
	while (vector[i] < 0)
	{
		sum += vector[i];
		vector[i] = 0;
		i++;
	}

	j = size - i;
	while (vector[i] + sum / j < 0)
	{
		sum += vector[i];
		vector[i] = 0;
		j--; i++;
	}

	for (i; i < size; i++)
		vector[i] += sum / j;

	double *temp = new double[size];
	for (int i = 0; i < size; temp[i] = vector[i++]);
	for (int i = 0; i < size; vector[order[i]] = temp[i++]);

	if (order != NULL) delete[]order;
	if (temp != NULL) delete[]temp;
}
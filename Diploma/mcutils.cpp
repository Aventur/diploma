#include "mcutils.h"

istringstream * get_string_stream(istream &is, const char * alphabet)
{
	is.seekg(0, ios::end);						// ������� ������������������
	int g = is.tellg();							// �� ����� � �����������
	is.seekg(0, ios::beg);						// ������
	char *file = new char[g];
	is.read(file, g);

	int new_g = 0;								// ������� ���������� ���������� ��������
	for (int i = 0; i < g; i++)
		if (strchr(alphabet, file[i]) != NULL)
			new_g++;
	char *new_file = new char[new_g];			// ����� ������, ��� ���������� ��������

	new_g = 0;
	for (int i = 0; i < g; i++)					// ������� ���������� ��������
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

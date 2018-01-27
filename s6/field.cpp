#include "field.h"
#include <cstdarg>
#include <iostream>
#include <cmath>

using namespace std;
using namespace poisson;

field::field(int D, int col): D(D), col(col)	
{
	int max = (int) pow((double)col, (double)D);
	campo = new double[max]();
}

field::~field()
{
	delete[] campo;
}

double field::getValue(int x0, ...)
{
	va_list va;
	va_start(va, x0);

	int offset = getOffset(x0, va);

	va_end(va);

	return campo[offset];
}

double field::getValue(int* coords)
{
	return campo[getOffset(coords)];	
}

void field::setValue(double v, int x0, ...)
{
	va_list va;
	va_start(va, x0);

	int offset = getOffset(x0, va);

	va_end(va);
	campo[offset] = v;
}

void field::setValue(double v, int* coords)
{
	campo[getOffset(coords)] = v;
}

int field::getOffset(int x0, va_list& va)
{
	int offset = x0;

	for(int i = 1; i<D; i++)
	{
		offset += (va_arg(va, int) * pow(this->col, i));
	}

	return offset;
}

int field::getOffset(int* coords)
{
	int offset = coords[0];
	for(int i = 1; i<D; i++)
	{
		offset += (coords[i] * pow(this->col, i));
	}

	return offset;

}

void field::fill(ifstream& in)
{
	int *coords = new int[D]();

	while(!in.eof())
	{
		for(int i = 0; i<D; i++)
		{
			in >> coords[i];
		}

		int offset = getOffset(coords);

		in >> campo[offset];
	}

	delete[] coords;
}

void field::print(ofstream& out)
{
	int *coords = new int[this->D]();
	int max = (int) pow((double)this->col, (double)this->D);

	for(int i = 0; i<max; i++)
	{

		int j = 0;
		for(j = 0; j<D; j++)
		{
			out << coords[j] << "\t";
			//cout << coords[j] << "\t";
		}

		int of = getOffset(coords);
		out << campo[of] << endl;
		//cout << campo[of] << endl;

		j = 0;
		while(j<D)
		{
			if(coords[j] < (this->col-1))
			{
				coords[j]++;
				break;
			}
			else
			{
				coords[j] = 0;
			}
			j++;
		}
	}
}



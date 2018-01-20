#pragma once

#include <cstdarg>
#include <fstream>

using namespace std;

namespace poisson
{
	class field
	{
	private:
		double* campo;
		int D;
		int col;

		int getOffset(int x0, va_list& va);
		int getOffset(int* coords);

	public:
		field(int D, int col);

		~field();

		double getValue(int x0, ...);
		double getValue(int* coords);
		void setValue(double v, int x0, ...);
		void setValue(double v, int* coords);
		void fill(ifstream& in);
		void print(ofstream& in);
		//void setCharge()
	};
}
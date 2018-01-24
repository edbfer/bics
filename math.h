#pragma once

#include "complex.h"
#include "matriz.h"

#include <fstream>

using namespace std;

namespace math
{
	double h;

	template <typename T>
	complex<T> simpsonPond(int i, int j, complex<T> vij)
	{
		/*double multi;
		if(	(i == 0 && j == 0) ||
				(i == 0 && j == 127) ||
				(i == 127 && j == 0) ||
				(i == 127 && j == 127))
		{
			multi = 1;
		} else if((i == 0 || j == 0) && (i != 127 && j != 127))
		{
			if(i % 2 == 0)
				multi = 2;
			else
				multi = 4;
		} else
		{
			if(i % 2 == 0)
			{
				if(j % 2 == 1)
					multi = 8;
				else
					multi = 4;
			}
			else
			{
				if(j % 2 == 1)
					multi = 16;
				else
					multi = 8;
			}
		}*/
		double multi = 0;
		if(	(i == 0 && j == 0) ||
				(i == 0 && j == 127) ||
				(i == 127 && j == 0) ||
				(i == 127 && j == 127))
		{
			multi = 1;
		}
		else if((i == 0 || i == 127) && (j != 0 && j != 127))
		{
			if(j % 2 == 0)
				multi = 2;
			else
				multi = 4;
		} else if((j == 0 || j == 127) && (i != 0 && i != 127))
		{
			if(i % 2 == 0)
				multi = 2;
			else
				multi = 4;
		}else
		{
			if(i % 2 == 0)
			{
				if(j % 2 == 0)
					multi = 4;
				else
					multi = 8;
			}
			else
			{
				if(j % 2 == 0)
					multi = 8;
				else
					multi = 16;
			}
		}

		return vij * multi;
	}

	template <typename T>
	complex<T> conjugate(int i, int j, complex<T> vij)
	{
		return vij * (~vij);
	}

	template <typename T>
	complex<T> simpson2d(matriz<T> f)
	{
		ofstream out("fixeru.txt");
		matriz<T> pond = f.execute(conjugate<T>);
		pond = pond.execute(simpsonPond<T>);
		complex<T> v = 0.;
		for (int i = 0; i<f.n; i++)
		{
			for (int j = 0; j<f.m; j++)
			{
				v = v + pond(i, j);
			}
		}

		v = v * ((1./9.)*h*h);
		return v;
	}

	template <typename T>
	matriz<T> dx(matriz<T> f)
	{
		matriz<T> res(f.n, f.m);
		matriz<float> d = matriz<float>::tridiagonal(-1, 0, 1, 128);
		//f = transpose(f);
		res = d * f;
		res = res * (complex<T>)(0.5*h);
		return res;
	}

	template <typename T>
	matriz<T> dy(matriz<T>& f)
	{
		matriz<T> res(f.n, f.m);
		matriz<float> d = matriz<float>::tridiagonal(-1, 0, 1, 128);

		res = d * f;
		res = res * (complex<T>)(1/(2*h));
		return res;
	}

	template <typename T>
	matriz<T> d2x(matriz<T>& f)
	{
		matriz<T> res(f.n, f.m);
		matriz<float> d2 = matriz<float>::tridiagonal(1, -2, 1, 128);
		res = d2 * f;
		res = res * (1/(h*h));
		return res;
	}

	template <typename T>
	matriz<T> d2y(matriz<T>& f)
	{
		matriz<T> res = (f.n, f.m);
		matriz<float> d2 = matriz<float>::tridiagonal(1, -2, 1, 128);
		res = transpose(f);
		res = d2 * res;
		res = res * (1/(h*h));
		return res;
	}
}

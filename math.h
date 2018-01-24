#pragma once

#include "complex.h"
#include "matriz.h"

#include <fstream>

using namespace std;

namespace math
{
	double h;

	template <typename T>
	complex<T> simpson2d(matriz<T> f)
	{
		ofstream out("fixeru.txt");
		matriz<T> pond(128,128);

		for(int i=0;i<128;i++)
		{
			if(i==0 || i==127)
			{
				pond(i,0)=(complex<T>)1;
				pond(i,127)=(complex<T>)1;
			}

			if(i!=0 && i!=127 && i%2==0)
			{
				pond(i,0)=(complex<T>)2;
				pond(i,127)=(complex<T>)2;
			}else
			{
				pond(i,0)=(complex<T>)4;
				pond(i,127)=(complex<T>)4;
			}
		}

		for(int j=1;j<127;j++)
		{
			for(int i=0;i<128;i++)
			{
				if(j%2==0)
				{
					pond(i,j)=((complex<T>)2)*pond(0, j);
				} else
				{
					pond(i,j)=((complex<T>)4)*pond(0, j);
				}
			}
		}

		out << pond;

		return 0;
	}

	template <typename T>
	matriz<T> dx(matriz<T> f)
	{
		matriz<T> res(f.n, f.m);
		matriz<float> d = matriz<float>::tridiagonal(-1, 0, 1, 128);
		//f = transpose(f);
		res = d * f;
		res = res * (complex<T>)(0.5);
		return res;
	}

	template <typename T>
	matriz<T> dy(matriz<T>& f)
	{
		matriz<T> res(f.n, f.m);
		matriz<float> d = matriz<float>::tridiagonal(-1, 0, 1, 128);

		res = d * f;
		//res = res * (complex<T>)(1/(0.5*h));
		return res;
	}
}

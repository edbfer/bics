#pragma once

#include "complex.h"
#include "field.h"
#include "matriz.h"

#include <fstream>

using namespace std;

namespace math
{
template <typename T>
complex<T> simpson2d(field<T> f)
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
					pond(i,j)=((complex<T>)2)*pond(i,1);
				} else
				{
					pond(i,j)=((complex<T>)4)*pond(i,1);
				}
			}
		}

		out << pond;

		return 0;
	}
}

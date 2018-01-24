#include "complex.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace math;

int main(int argc, char const *argv[])
{
	double gamma, omega, G;

	matriz<float> lattice(512, 512);

	ifstream f("params.txt");
	f >> gamma >> omega >> G;
	f.close();

	cout << "Condição inicial: " << endl;
	float cond;
	cin >> cond;

	cout << "Tempo máximo: " << endl;
	float tmax;
	cin >> tmax;

	cout << "dt: " << endl;
	float dt;
	cin >> dt;

	/*lattice.fill(cond);

	ifstream b("cond.txt");
	lattice.fill(b);
	b.close();

	ofstream init("initial.txt");
	lattice.printCoord(init);
	init.close();	*/
	math::h = 40.0/512.0;

	double x, y;
	for (int i = 0; i<512; i++)
	{
		for (int j = 0; j<512; j++)
		{
			x = -20 + h*i;
			y = -20 + h*j;
			lattice(i, j) = cos(sqrt((x*x + y*y)));
		}
	}

	ofstream init("initial.txt");
	lattice.printCoord(init);
	init.close();

	//complex<float> norm = simpson2d<float>(lattice);
	lattice = dx<float>(lattice);
	ofstream ddx("dx.txt");
	lattice.printCoord(ddx);
	ddx.close();

	lattice = dy<float>(lattice);
	ofstream res("dy.txt");
	lattice.printCoord(res);

	return 0;
}
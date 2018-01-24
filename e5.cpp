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

	matriz<float> lattice(128, 128);

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

	/*lattice.fill(cond);2acf2de1-444c-4b39-b707-323917289e24

	ifstream b("cond.txt");
	lattice.fill(b);
	b.close();

	ofstream init("initial.txt");
	lattice.printCoord(init);
	init.close();	*/
	math::h = 20.0/129.0;

	double x, y;
	for (int i = 0; i<128; i++)
	{
		for (int j = 0; j<128; j++)
		{
			x = -10 + h*i;
			y = -10 + h*j;
			lattice(i, j) = 1/x;
		}
	}

	complex<float> norm = simpson2d<float>(lattice);

	/*ofstream init("initial.txt");
	init << lattice;*/

	//complex<float> norm = simpson2d<float>(lattice);

	cout << "Integral: " << norm << endl;

	/*lattice = dx<float>(lattice);
	ofstream ddx("dx.txt");
	lattice.printCoord(ddx);
	ddx.close();

	lattice = dy<float>(lattice);
	ofstream res("dy.txt");
	lattice.printCoord(res);*/

	return 0;
}

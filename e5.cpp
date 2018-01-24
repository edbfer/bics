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

	lattice.fill(cond);

	ifstream b("cond.txt");
	lattice.fill(b);
	b.close();

	ofstream init("initial.txt");
	lattice.printCoord(init);
	init.close();	

	math::h = 20.0/128.0;

	//complex<float> norm = simpson2d<float>(lattice);
	lattice = dx<float>(lattice);
	//lattice = dy<float>(lattice);

	ofstream res("result.txt");
	lattice.printCoord(res);

	return 0;
}
#include "complex.h"
#include "field.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace math;

int main(int argc, char const *argv[])
{
	double gamma, omega, G;

	field<float> lattice(2, 128);

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

	complex<float> norm = simpson2d<float>(lattice);

	matriz<float> tri = matriz<float>::tridiagonal(1, 0, 1, 128);
	ofstream print("outma.txt");
	print << tri;

	ofstream res("result.txt");
	lattice.print(res);

	return 0;
}
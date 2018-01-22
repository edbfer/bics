#include "complex.h"
#include "field.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

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

	lattice.fill(cond);

	ifstream b("cond.txt");
	lattice.fill(b);
	b.close();

	complex<float> norm = simpson2d(lattice);

	ofstream res("result.txt");
	lattice.print(res);

	return 0xdeadbeef;
}
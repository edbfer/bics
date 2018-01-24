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

	complex<float> i(0, 1);

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

	init.close();
	init.close();	*/
	math::h = 20.0/129.0;
	 derivuxa= dx<float>(lattice);
	 derivuya= dy<float>(lattice);
	//lattice = dy<float>(lattice);


	cout << "Integral: " << norm << endl;

	/*lattice = dx<float>(lattice);
	ofstream ddx("dx.txt");
	lattice.printCoord(ddx);
	ddx.close();

	lattice = dy<float>(lattice);
	ofstream res("dy.txt");
	lattice.printCoord(res);*/

	matriz<float> lapx= d2x<float> (lattice);
	matriz<float> lapy= d2y <float>(lattice);
	matriz<float> lap = lapx + lapy;

	matriz<float> laplacianu(128,128)= (-0.5*lap);
	matriz<float> paressela2(128,128)= ((x*x + y*y)/2)*lattice;
	matriz<float> paressela3(128,128)= G*norm*norm*lattice;

	matriz<float> parssialex(128,128)=y*dx;
	matriz<float> parssialey(128,128)=x*dy;

	matriz<float> ultimeparcele(128,128)= i*omega*(parssialey-parssialex);


	return 0;
}

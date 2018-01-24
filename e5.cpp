#include "complex.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;
using namespace math;

template <typename T>
complex<T> x2y2(int i, int j, complex<T> vij)
{
	double x = math::x0 + i*math::h;
	double y = math::y0 + j*math::h;
	return ((0.5)*(x*x + y*y)) * vij;
}

template <typename T>
complex<T> valx(int i, int j, complex<T> vij)
{
	double x = math::x0 + i*math::h;
	return x * vij;
}

template <typename T>
complex<T> valy(int i, int j, complex<T> vij)
{
	double y = math::y0 + i*math::h;
	return y * vij;
}

template <typename T>
matriz<T> F(double G, double gamma, double omega, matriz<T>& lattice)
{
		complex<float> i(0, 1);

		matriz<float> lapx= d2x<float> (lattice);
		matriz<float> lapy= d2y <float>(lattice);
		matriz<float> lap = lapx + lapy;

		matriz<float> laplacianu = lap * (-0.5);

		matriz<float> paressela2 = lattice.execute(x2y2);

		complex<float> norm = simpson2d<float>(lattice);
		matriz<float> paressela3 = lattice*(G*norm*norm);

		matriz<float> ddx = math::dx(lattice);
		matriz<float> ddy = math::dy(lattice);
		ddx = ddx.execute(valy<float>);
		ddy = ddy.execute(valx<float>);

		matriz<float> ultimeparcele = ddx + ddy;
		ultimeparcele = ultimeparcele * (i * omega);

		matriz<float> total = laplacianu + paressela2 + paressela3 + ultimeparcele;
		total = total * (1/(i - gamma));

		return total;
}


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

	lattice.fill(1.);
	ifstream c("cond.txt");
	lattice.fill(c);
	c.close();

	cout << "Tempo máximo: " << endl;
	float tmax;
	cin >> tmax;

	cout << "dt: " << endl;
	float dt;
	cin >> dt;

	math::h = 20.0/129.0;

	complex<float> norm = simpson2d<float>(lattice);
	lattice = lattice * (1/sqrt((float) norm));

	int i = 0;

	for(double t = 0; t<tmax; t += dt)
	{
		matriz<float> flat = F<float>(G, gamma, omega, lattice);

		matriz<float> psi1 = lattice + (flat * dt);
		matriz<float> fpsi1 = F<float>(G, gamma, omega, psi1);

		matriz<float> psi2 = (lattice*(3./4.)) + (psi1*(1./2.)) + (fpsi1*(dt*(1./4.)));
		matriz<float> fpsi2 = F<float>(G, gamma, omega, psi2);

		matriz<float> psitotal = (lattice*(1./3.)) + (psi2*(2./3.)) + (fpsi2*(dt*(2./3.)));

		complex<float> norm = simpson2d<float>(psitotal);
		psitotal = psitotal * (1/sqrt((float) norm));

		stringstream file;
		file << "dados/t" << t << ".txt";
		ofstream out(file.str());

		psitotal.printCoord(out);

		out.close();

		lattice	= psitotal;

		cout << "Iteração " << ++i << endl;
	}
	return 0;
}

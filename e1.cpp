#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <omp.h>
#include <time.h>
#include <sstream>

//#include "cuda.cuh"

#include "matriz.h"
#include "complex.h"

using namespace std;

int main(int argc, char const *argv[])
{
	//matriz m1(3, 3);
	//cin >> m1;
	srand(time((time_t*) 0));

	/*matriz<complex> m1 = matriz<complex>(3, 3);
	cin >> m1;
	matriz<complex> m2 = m1.eigenpairs(m1, 100);

	ofstream o("out.txt");
	o << fixed << setw(10) << setprecision(8);
	o << m2 << endl;

	system("pause");

	return 0;*/

	int dim;
	double l;
	cout << "Inserir Dimensão Matricial: " << endl;
	cin >> dim;
	cout << "Inserir intervalo: " << endl;
	cin >> l;

	double h = (2*l)/(dim + 1);

	cout << "H: " << h << endl;

	matriz K = matriz::id(dim);
	matriz V = matriz::id(dim);

	#pragma omp parallel for
	for(int i = 0; i<dim-1; i++)
	{
		K(i, i) = -2;
		K(i + 1, i) = 1;
		K(i, i + 1) = 1;
	}
	K(dim-1, dim-1) = -2;


	double x0 = -l;
	#pragma omp parallel for
	for(int i = 0; i<dim; i++)
	{
		double x = x0 + i*h;
		V(i, i) = (1./2.)*x*x;
	}

	K = K * (complex) (-1./2.);

	matriz H = K + V;

	ofstream oH("schrodinger/h.txt");
	oH << H << endl;
	oH.close();

	clock_t t = clock();

	ofstream oRes("schrodinger/res.txt");
	matriz resultados = matriz::eigenpairs(H, 100);

	t = clock() - t;
	oRes << resultados << endl;
	oRes.close();

	cout << "Time elapsed: " << ((double)t/CLOCKS_PER_SEC) << endl;

	#pragma omp parallel for
	//merda pa-
	for(int i = 0; i<dim; i++)
	{
		ostringstream os;
		os << "schrodinger/psi_n" << i << ".txt";
		ofstream out(os.str().c_str());
		matriz v = extract(resultados, 0, dim-1, i, i);
		const double x0 = -l;
		for(int j = 0; j<dim; j++)
		{
			double x = x0 + j*h;
			out << x << "\t" << v(j, 0) << endl;
		}

	}

	return 0;
	
	/*
	matriz m2(dim, dim);
	char res;
	cout << "Auto preencher? " << endl;
	cin >> res;

	if (res == 'Y' || res == 'y')
	{
		m2 = matriz::hermitian(dim);
	}
	else
	{
		cin >> m2;
	}

	cout << fixed << setw(10) << setprecision(8) << endl;

	ofstream out("det.txt");
	double x = 0.0f;
	matriz i = matriz::id(dim);
	#pragma omp parallel for
	for(; x<4.0f; x+=1e-5 )
	{
		matriz m3 = m2-(i*x);
		double det = determinant(m3);
		out << x << "\t" << det << endl;
	}

	matriz m4 = matriz::eigenpairs(m2, 100);
	cout << m4 << endl;

	return 0;
	/* //antigo código*/
}


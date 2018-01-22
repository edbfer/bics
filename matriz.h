#ifndef _MATRIZ_H
#define _MATRIZ_H 1

#include <cmath>
#include <iostream>
#include "complex.h"
#include <vector>

using namespace std;

template <typename T>
class matriz
{

	struct swapOut
	{
		complex pivot;
		int swapped;
		int l1;
		int l2;
	};

	struct gaussInfo
	{
		int ispermutation;	
		int p1;
		int p2;
		complex v;
	};

public:
	//complex ** mat;
	complex * mat;
	int n, m;
	mutable complex det;
	mutable int dirty;
	mutable int sigma;

	matriz(const int n = 3, const int m = 3);
	matriz(const matriz& m1);
	~matriz();

	static matriz id(int n);
	friend matriz permutation(int n, int l1, int l2);
	static matriz random(int n, int m);
	static matriz hermitian(int n);

	inline complex& operator()(const int& n, const int& m)
	{
		return mat[n * this->m + m];
	}

	complex* operator[](const size_t i) const;
	matriz& operator=(const matriz& m1);

	matriz operator-();
	matriz operator++(int);
	matriz operator~();

	friend matriz transpose(matriz& m1);
	static matriz inverte(matriz& m1);
	friend matriz extend(matriz& m1, matriz& m2);
	friend matriz extract(matriz& m1, int x0, int y0, int x1, int y1);

	void swapLine(int l1, int l2);
	matriz swapColForVector(int col, matriz& vec);

	friend matriz operator^(const matriz& m1, const int e);
	friend matriz operator+(const matriz& m1, const matriz& m2);
	//friend matriz sum_cuda(matriz& m1, matriz& m2);
	friend matriz operator-(const matriz& m1, const matriz& m2);
	friend matriz operator*(const matriz& m1, const matriz& m2);
	friend matriz operator~(const matriz& m1);
	friend matriz multiply_cuda(matriz& m1, matriz& m2);
	friend matriz operator*(const matriz& m1, complex v);

	friend complex dot(matriz& m1, matriz& m2);

	friend ostream& operator<<(ostream& out, matriz& m1);
	friend istream& operator>> (istream& in, matriz& m1);

	friend complex determinant(matriz& m1);
	swapOut swapLineForPivot(int col, int inv);
	friend complex normInf(matriz& m1);

	matriz gauss(int inv, int c1, int c2);
	friend matriz gramschmidt(matriz& m1);

	matriz LUDecomposition(int inv);
	matriz QRDecomposition();

	static matriz eigenpairs(matriz& m1, int nmax);


};

#endif //_MATRIX_H
#ifndef _MATRIZ_H
#define _MATRIZ_H 1

#include <cmath>
#include <iostream>
#include "complex.h"
#include <vector>
#include <cstring>

using namespace std;

template <typename T>
class matriz
{

	struct swapOut
	{
		complex<T> pivot;
		int swapped;
		int l1;
		int l2;
	};

	struct gaussInfo
	{
		int ispermutation;	
		int p1;
		int p2;
		complex<T> v;
	};

public:
	//complex ** mat;
	complex<T> * mat;
	int n, m;
	mutable complex<T> det;
	mutable int dirty;
	mutable int sigma;

	matriz(const int n = 3, const int m = 3)
	{
		mat = new complex<T>[n*m];
		dirty = 1;
		det = (T)0.;
	}
	matriz(const matriz& m1)
	{
		mat = new complex<T>[n*m];
		memcpy(mat, m1.mat, sizeof(complex<T>)*n*m);
		det = m1.det;
		dirty = m1.dirty;
	}
	~matriz()
	{
		delete[] mat;
	}

	static matriz id(int n)
	{
		matriz res(n, n);
		for(int i = 0; i<n; i++)
		{
			res(i, i) = ((complex<T>) 1.);
		}
		return res;
	}
	friend matriz permutation(int n, int l1, int l2)
	{
		matriz res = matriz::id(n);
		res.swapLine(l1, l2);
		return res;
	}
	static matriz random(int n, int m)
	{
		matriz res(n, m);
		for(int i = 0; i<n; i++)
		{
			for(int j = 0; j<m; j++)
			{
				res(i, j) = (complex<T>) ((double)rand()/(double) RAND_MAX)*10;
			}
		}
		return res;
	}
	static matriz hermitian(int n)
	{
		matriz res(n, n);
		for(int i = 0; i<n; i++){
			for(int j = 0; j<i+1; j++)
			{
				res(i, j) = (complex<T>) ((double)rand()/(double) RAND_MAX)*10;
				res(j, i) = res(i, j);
			}
		}
		return res;
	}
	static matriz tridiagonal(complex<T> d0, complex<T> dp, complex<T> d1, int x)
	{
		matriz res(x, x);
		for (int i = 0; i<x-1; i++)
		{
			res(i, i) = dp;
			res(i, i + 1) = d1;
			res(i + 1, i) = d0;
		}
		res(x-1, x-1) = dp;
		return res;
	}

	inline complex<T>& operator()(const int& n, const int& m)
	{
		return mat[n * this->m + m];
	}

	complex<T>* operator[](const size_t i) const
	{
		return (mat + i*m);
	}
	matriz& operator=(const matriz& m1)
	{
		delete[] mat;

		n = m1.n;
		m = m1.m;

		mat = new complex<T>[n*m];
		memcpy(mat, m1.mat, sizeof(complex<T>)*n*m);

		dirty = m1.dirty;
		det = m1.det;

		return *this;
	}

	matriz operator-()
	{
		for(int i = 0; i<n; i++)
		{
			for(int j = 0; j<m; j++)
			{
				(*this)(i, j) = -(*this)(i, j);
			}
		}
		return *this;
	}
	matriz operator++(int)
	{
		for(int i = 0; i<n; i++)
		{
			for(int j = 0; j<m; j++)
			{
				(*this)(i, j) = (*this)(i, j) + (complex<T>) 1.;
			}
		}
		return *this;
	}

	friend matriz transpose(matriz& m1)
	{
		matriz r(m1.m, m1.n);
		for(int i = 0; i<m1.n; i++)
		{
			for(int j = 0; j<m1.m; j++)
			{
				r(j, i) = m1(i, j);
			}
		}
		return r;
	}
	static matriz inverte(matriz& m1)
	{
		complex<T> det = determinant(m1);
		if(det == 0.)
		{
			cout << "Singular" << endl;
			exit(-1);
		}
		cout << "Determinante: " << det  << "\n\n\n\n" << endl;


		matriz i = matriz::id(m1.m);
		matriz t = extend(m1, i);
	//cout << t << endl;

		matriz res = t.gauss(0, 0, m1.n);
	//cout << res << endl;

		res = res.gauss(1, 0, m1.n-1);
	//cout << res << endl;


		for(int i = 0; i<m1.n; i++)
		{
			complex<T> v = res(i, i);
			for(int j = 0; j<res.m; j++)
			{
				res(i, j) = res(i, j) / v;
			}
		}

		res = extract(res, 0, m1.n-1, m1.m, t.m-1);

		return res;
	}
	friend matriz extend(matriz& m1, matriz& m2)
	{
		matriz t(m1.n, m1.m+m2.m);
		for(int i = 0; i<t.n; i++)
		{
			for(int j = 0; j<m1.m; j++)
			{
				t(i, j) = m1(i, j);
			}

			for(int j = 0; j<m2.m; j++)
			{
				t(i, j+m1.m) = m2(i, j);
			}
		}

		return t;
	}
	friend matriz extract(matriz& m1, int x0, int y0, int x1, int y1)
	{
		matriz r(x1 - x0 + 1, y1 - y0 + 1);
		for (int i = x0; i <=x1; i++)
		{
			for (int j = y0; j <= y1; j++)
			{
				r(i - x0, j - y0) = m1(i, j);
			}
		}
		return r;
	}

	void swapLine(int l1, int l2)
	{
		complex<T>* r = new complex<T>[m];
		memcpy(r, (*this)[l1], sizeof(complex<T>)*m);
		memcpy((*this)[l1], (*this)[l2], sizeof(complex<T>)*m);
		memcpy((*this)[l2], r, sizeof(complex<T>)*m);
	}
	matriz swapColForVector(int col, matriz& vec)
	{
		for (int i = 0; i < n; i++)
		{
			(*this)(i, col) = vec(i, 0);
		}
		return *this;
	}

	friend matriz operator^(const matriz& m1, const int e)
	{
		matriz r = m1;
		for(int i = 0; i<e; i++)
		{
			r = r * m1;
		}
		return r;
	}
	friend matriz operator+(const matriz& m1, const matriz& m2)
	{
		matriz res(m1.n, m1.m);
		for(int i = 0; i<m1.n; i++)
		{
			for(int j = 0; j<m1.m; j++)
			{
				res(i, j) = m1(i, j) + m2(i, j);
			}
		}
		return res;
	}
	//friend matriz sum_cuda(matriz& m1, matriz& m2);
	friend matriz operator-(const matriz& m1, const matriz& m2)
	{
		matriz res(m1.n, m1.m);
		for(int i = 0; i<m1.n; i++)
		{
			for(int j = 0; j<m1.m; j++)
			{
				res(i, j) = m1(i, j) - m2(i, j);
			}
		}
		return res;
	}
	friend matriz operator*(const matriz& m1, const matriz& m2)
	{
		matriz newm(m1.n, m2.m);
		for (int i = 0; i<newm.m; i++)
		{
			for(int j = 0; j<newm.n; j++)
			{
				complex<T> v = 0.;
				for(int k = 0; k<m1.m; k++)
				{
					complex<T> v1 = m2(k, i);
					complex<T> v2 = m1(j, k);
					v = v + v1 * v2;
				}
				newm(j, i) = v;
			}
		}

		return newm;
	}
	matriz operator~()
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				~(*this)(i, j);
			}
		}
		return *this;
	}
	//friend matriz multiply_cuda(matriz& m1, matriz& m2);
	friend matriz operator*(const matriz& m1, complex<T> v)
	{
		matriz res = m1;
		int i = 0, j = 0;
		for(; i<m1.n; i++)
		{
			for(j = 0; j<m1.m; j++)
			{
				res(i, j) = res(i, j) * v;
			}
		}

		return res;
	}

	friend complex<T> dot(matriz& m1, matriz& m2)
	{
		matriz temp = transpose(m1)*(~m2);
		return temp(0, 0);
	}

	friend ostream& operator<<(ostream& out, matriz& m1)
	{
		for (int i = 0; i<m1.n; i++)
		{
			for(int j = 0; j<m1.m; j++)
			{
				out << m1(i, j) << " ";
			}

			out << endl;
		}
		return out;
	}
	friend istream& operator>> (istream& in, matriz& m1)
	{
		for (int i = 0; i<m1.n; i++)	
		{
			for(int j = 0; j<m1.m; j++)
			{
				in >> m1(i, j);
			}
		}
		return in;
	}

	friend complex<T> determinant(matriz& m1)
	{
		if(!m1.dirty)
			return m1.det;
		else
		{
			m1.dirty = 0;
			m1.gauss(0, 0, m1.n - 1);
			return m1.det;
		}
	}
	swapOut swapLineForPivot(int col, int inv)
	{
		int firstZero = -1;
		int notZero = -1;
		int inc = (inv)?-1:+1;
		int i;
		for (i = col; (inv)?(i>=col):(i<m); i += inc)
		{
			if((*this)(i, col) != 0.)
				break;

		}

		swapOut out = {(*this)(i, col), (i == col)?0:1, i, col};

		swapLine(i, col);

		return out;
	}
	friend complex<T> normInf(matriz& m1)
	{
		vector<complex<T>> v;
		complex<T> s = 0.;
		for(int i = 0; i<m1.n; i++)
		{
			s = 0.;
			for(int j = 0; j<m1.m; j++)
			{
				s = s + m1(i, j);
			}
			v.push_back(s);
		}

		sort(v.begin(), v.end());

		return v[v.size()-1];
	}

	matriz gauss(int inv, int c1, int c2)
	{
		matriz m1 = *this;
		int inc = (inv)?-1:1;
		sigma = 1.;
		for(int col = (inv)?c2:c1; (inv)?col>=c1:col<c2; col += inc)
		{
			swapOut inf = m1.swapLineForPivot(col, inv);
			sigma = (inf.swapped)?sigma*-1:sigma;

			for(int row = (inv)?col-1:col+1; (inv)?row>=c1:row<c2; row += inc)
			{
				complex<T> v1;
				if((v1 = m1(row, col)) == 0)
					continue;

				complex<T> val = v1 / inf.pivot;
				for(int k = 0; k<m; k++)
				{
					m1(row, k) = m1(row, k) - (val)*m1(col, k);
				}

			}
		}

		this->det = 1.;
		for(int i = 0; i<n; i++)
		{
			this->det = det * (*this)(i, i);
		}

		this->det = det * sigma;
		m1.det = det;

		return m1;
	}
	friend matriz gramschmidt(matriz& m1)
	{
		matriz res = m1;

		for(int i = 0; i<m1.m; i++)
		{
			matriz v = extract(m1, 0, m1.n-1, i, i);
			for(int j = 0; j<i; j++)
			{
				matriz u = extract(res, 0, res.n-1, j, j);
			/*complex vu = dot(u, v);
			complex mod = dot(u, u);*/
				matriz temp = u*((complex<T>) (dot(u, v) / dot(u, u)));
				v = v - temp;
			}
			complex<T> mod = (complex<T>) sqrt((double) dot(v, v));
			v = v * (complex<T>) (1 / mod);
			res.swapColForVector(i, v);
		}
		return res;
	}

	matriz LUDecomposition(int inv)
	{
		vector<gaussInfo> v;
		complex<T> sigma = 1.;
	//cout << r << endl;
		for(int col = (inv)?m-1:0; (inv)?col>=0:col<m; (inv)?col--:col++)
		{
			swapOut inf = swapLineForPivot(col, inv);
			gaussInfo i = {1, inf.l1, inf.l2, 0.};
			v.push_back(i);
			sigma = (inf.swapped == 1)?sigma*-1:sigma*1;

			for(int row = (inv)?col-1:col+1; (inv)?row>=0:row<n; (inv)?row--:row++)
			{
				complex<T> v1 = (*this)(row, col);	

				if(v1 == 0.)
					continue;

				for(int k = col; (inv)?k>=0:k<m; (inv)?k--:k++)
				{
					complex<T> val = (*this)(row, k);
					complex<T> l1v = (*this)(col, k);
					(*this)(row, k) = val - (v1/inf.pivot)*l1v;
				}
				gaussInfo e = {0, row, col, -(v1/inf.pivot)};
				v.push_back(e);
			}

		}

		matriz res = matriz::id(n);
		for (int i = v.size()-1; i>=0; i--)
		{
			if (v[i].ispermutation == 1)
			{
				matriz p = permutation(n, v[i].p1, v[i].p2);
				res = res * p;
			}
			else
			{
				matriz e = matriz::id(n);
				e(v[i].p1, v[i].p2) = v[i].v;
				res = res * e;
			}
		}

		det = 1.;
		for(int i = 0; i<n; i++)
		{
			det = det * (*this)(i, i);
		}

		det = det * sigma;
		return res;

	}
	matriz QRDecomposition()
	{
		matriz A = *this;
		*this = gramschmidt(*this);
		matriz qinv = transpose(*this);
		matriz R = qinv * A;
		return R;
	}

	static matriz eigenpairs(matriz& m1, int nmax)
	{
		matriz temp = m1;
		matriz v0 = random(m1.n, m1.m);
		v0 = gramschmidt(v0);
		matriz v1 = matriz::id(m1.n);

		for (int i = 0; i < nmax; i++)
		{
			v1 = multiply_cuda(temp, v0);
			v0 = gramschmidt(v1);
			cout << "Iteração: " << i << " de " << nmax << endl;
		}

		matriz eigenvalues(m1.n, 1);
	//rayleigh quotien
		for (int i = 0; i < m1.m; i++)
		{
			matriz x = extract(v0, 0, v0.n - 1, i, i);
			complex<T> t = (complex<T>) ((transpose(x)*temp)*x)(0, 0);
			complex<T> d = (complex<T>) dot(x, x);
			eigenvalues(i, 0) = t/d;
		}

		v0 = extend(v0, eigenvalues);
		return v0;
	}
};

#endif //_MATRIX_H
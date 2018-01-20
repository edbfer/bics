#pragma once

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace std;

template <typename T>
class complex
{
	T real;
	T im;

public:
	__host__ __device__ complex() : real(0.), im(0.) {}
	__host__ __device__ complex(int real) : real(real), im(0.) {}
	__host__ __device__ complex(T real) : real(real), im(0.) {}
	__host__ __device__ complex(T real, T im) : real(real), im(0.) {}
	
	
	__host__ __device__ friend complex operator+(complex& c1, complex& c2)
	{
	  complex res(0, 0);
	  res.real = c1.real + c2.real;
	  res.im = c1.im + c2.im;
	  return res;
	}
	
	__host__ __device__ friend complex operator-(complex& c1, complex& c2)
	{
	  complex res(0, 0);
	  res.real = c1.real - c2.real;
	  res.im = c1.im - c2.im;
	  return res;
	}

	__host__ __device__ complex operator-()
	{
	  this->real *= -1;
	  this->im *= -1;
	  return *this;
	}
	
	__host__ __device__ complex operator~()
	{
	  im = -im;
	  return *this;
	}
	
	__host__ __device__ friend complex operator*(complex& c1, complex& c2)
	  {
	    complex res(0, 0);
	    res.real = c1.real * c2.real - c1.im * c2.im;
	    res.im = c1.im * c2.real + c1.real * c2.im;
	    return res;
	  }
	
	__host__ __device__ friend complex operator*(T a, complex& c1);
	{
	  complex res(0, 0);
	  res.real = a * c1.real;
	  res.im = a * c1.im;
	  return res;
	}
	
	__host__ __device__ friend complex operator/(complex& c1, complex& c2);
	{
	  complex res(0, 0);
	  T m = c2.real * c2.real + c2.im * c2.im;
	  res.real = (c1.real * c2.real + c1.im * c2.im)/m;
	  res.im	= (c1.im * c2.real - c1.real * c2.im)/m;
	  return res;
	}
	
	__host__ __device__ friend complex operator/(complex& c1, T a);
	{
	  complex res(0, 0);
	  res.real = c1.real/a;
	  res.im = c1.im/a;
	  return res;
	}

	__host__ __device__ friend int operator==(complex& c1, complex& c2);
	{
	  if ((c1.real == c2.real) && (c1.im == c2.im))
	    return 1;
	  return 0;
	}

	__host__ __device__ operator T()
	{
	  return real;
	}

	__host__ __device__ friend ostream& operator<<(ostream& out, complex& c1);
	{
	  const char* sign = (c1.im > 0) ? "+" : ((c1.im < 0) ? "-" : "+");
	  out << c1.real << sign << fabs(c1.im) << "i";
	  return out;
	}
	
	__host__ __device__ friend istream& operator>>(istream& in, complex& c1);
	{
	  in >> c1.real >> c1.im;
	  return in;
	}


	T mod();
};

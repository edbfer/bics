#pragma once

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace std;

class complex
{
	double real;
	double im;
	static int format;

public:
	__host__ __device__ complex();
	__host__ __device__ complex(int real);
	__host__ __device__ complex(double real);
	__host__ __device__ complex(double real, double im);
	
	__host__ __device__ friend complex operator+(complex& c1, complex& c2);
	__host__ __device__ friend complex operator-(complex& c1, complex& c2);
	__host__ __device__ complex operator-();
	__host__ __device__ complex operator~();
	__host__ __device__ friend complex operator*(complex& c1, complex& c2);
	__host__ __device__ friend complex operator*(double a, complex& c1);
	__host__ __device__ friend complex operator/(complex& c1, complex& c2);
	__host__ __device__ friend complex operator/(complex& c1, double a);

	__host__ __device__ friend int operator==(complex& c1, complex& c2);

	__host__ __device__ operator double()
	{
		return real;
	}

	__host__ __device__ friend ostream& operator<<(ostream& out, complex& c1);
	__host__ __device__ friend istream& operator>>(istream& in, complex& c1);

	double mod();
};
#pragma once

#include <iostream>

using namespace std;

template <typename T>
class complex
{
  T real;
  T im;

 public:
 complex() : real(0.), im(0.) {}
 complex(int real) : real(real), im(0.) {}
 complex(T real) : real(real), im(0.) {}
 complex(T real, T im) : real(real), im(0.) {}

  complex operator-()
  {
    this->real *= -1;
    this->im *= -1;
    return *this;
  }
	
  complex operator~()
  {
    im = -im;
    return *this;
  }
	
	
  friend complex operator+(complex& c1, complex& c2)
  {
    complex res(0, 0);
    res.real = c1.real + c2.real;
    res.im = c1.im + c2.im;
    return res;
  }
	
  friend complex operator-(complex& c1, complex& c2)
  {
    complex res(0, 0);
    res.real = c1.real - c2.real;
    res.im = c1.im - c2.im;
    return res;
  }

  friend complex operator*(complex& c1, complex& c2)
    {
      complex res(0, 0);
      res.real = c1.real * c2.real - c1.im * c2.im;
      res.im = c1.im * c2.real + c1.real * c2.im;
      return res;
    }
    

  friend complex operator*(T a, complex& c1)
    {
      complex res(0, 0);
      res.real = a * c1.real;
      res.im = a * c1.im;
      return res;
    }
	
  friend complex operator/(complex& c1, complex& c2)
  {
    complex res(0, 0);
    T m = c2.real * c2.real + c2.im * c2.im;
    res.real = (c1.real * c2.real + c1.im * c2.im)/m;
    res.im	= (c1.im * c2.real - c1.real * c2.im)/m;
    return res;
  }
	
  friend complex operator/(complex& c1, T a)
  {
    complex res(0, 0);
    res.real = c1.real/a;
    res.im = c1.im/a;
    return res;
  }

     
  friend int operator==(complex& c1, complex& c2)
  {
    if ((c1.real == c2.real) && (c1.im == c2.im))
      return 1;
    return 0;
  }


  operator T()
  {
    return real;
  }

  friend ostream& operator<<(ostream& out, complex& c1)
  {
    const char* sign = (c1.im > 0) ? "+" : ((c1.im < 0) ? "-" : "+");
    out << c1.real << sign << fabs(c1.im) << "i";
    return out;
  }
	
  friend istream& operator>>(istream& in, complex& c1)
  {
    in >> c1.real >> c1.im;
    return in;
  }


  T mod();
};


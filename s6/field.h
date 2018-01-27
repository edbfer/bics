#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include "matriz.h"

using namespace std;

template <typename T>
class field : public matriz<T>{
public:
  field(int col): matriz<T>(col, col){}

  field(const field<T>& f): matriz<T>(f)
  {}

  int getOffset(int* coords)
  {
    int D = 2;
    int offset = coords[0];
    for(int i = 1; i<D; i++)
    {
      offset += (coords[i] * pow(this->col, i));
    }
    return offset;
  }

  void fill(ifstream& in)
  {
    int *coords = new int[2]();

    while(!in.eof())
    {
      for(int i = 0; i<2; i++)
      {
        in >> coords[i];
      }

      int offset = getOffset(coords);
      in >> mat[offset];
    }

    delete[] coords;
  }

  void fill(T val)
  {
    int max = n*m;
    for(int i = 0; i < max; i++)
      mat[i] = val;
  } 
};

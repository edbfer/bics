#pragma once

#include <cstdarg>
#include <fstream>

using namespace std;

template <typename T>
class field
{
 private:
  typename T* campo;
  int D;
  int col;

  int getOffset(int* coords);

 public:
  field(int D, int col);

  ~field();

	        
  T getValue(int* coords);
  void setValue(T v, int* coords);
  void fill(ifstream& in);
  void print(ofstream& in);
  //void setCharge()
};

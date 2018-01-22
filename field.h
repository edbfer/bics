#pragma once

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

template <typename T>
class field
{
 private:
  T* campo;
  int D;
  int col; 

 public:
  field(int D, int col): D(D), col(col)  
  {
    int max = (int) pow((double)col, (double)D);
    campo = new T[max]();
  }

  ~field()
  {
    delete[] campo;
  }

  T getValue(int* coords)
  {
    return campo[getOffset(coords)];	
  }

  void setValue(T v, int* coords)
  {
    campo[getOffset(coords)] = v;
  }

  int getOffset(int* coords)
  {
    int offset = coords[0];
    for(int i = 1; i<D; i++)
      {
	     offset += (coords[i] * pow(this->col, i));
      }

    return offset;
  }

  
  void fill(ifstream& in)
  {
    int *coords = new int[D]();

    while(!in.eof())
      {
	for(int i = 0; i<D; i++)
	  {
	    in >> coords[i];
	  }

	int offset = getOffset(coords);

	in >> campo[offset];
      }

    delete[] coords;
  }

  void fill(T val)
  {
    int max = (int) pow((double)col, (double)D);
    for(int i = 0; i < max; i++)
      campo[i] = val;
  }


  void print(ofstream& out)
  {
    int *coords = new int[this->D]();
    int max = (int) pow((double)this->col, (double)this->D);

    for(int i = 0; i<max; i++)
      {

	int j = 0;
	for(j = 0; j<D; j++)
	  {
	    out << coords[j] << "\t";
	  }

	int of = getOffset(coords);
	out << campo[of] << endl;

	j = 0;
	while(j<D)
	  {
	    if(coords[j] < (this->col-1))
	      {
		coords[j]++;
		break;
	      }
	    else
	      {
		coords[j] = 0;
	      }
	    j++;
	  }
      }
  }  
};

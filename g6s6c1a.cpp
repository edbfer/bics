#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

//#include "cuda.cuh"

#include "complex.h"

using namespace std;


int main ()
{
  double DECISAOAOAOAO;
  double c;
  
  complex <double> c1;
  complex <double> c2;
  complex <double> result;
  
  cout << "Número complexo 1" << endl;
  cin >> c1;
  cout << endl;


  cout << "Oq keres fzr migu?" << endl;
  cout << "1) Soma" << endl;
  cout << "2) Subtração" << endl;
  cout << "3) Multiplicação (entre complexos)" << endl;
  cout << "4) Multiplicação (por escalar)" << endl;
  cout << "5) Divisão (entre complexos)" << endl;
  cout << "6) Divisão (por escalar)" << endl;

  cin >> DECISAOAOAOAO;

  if(DECISAOAOAOAO!=1 && DECISAOAOAOAO!=2 && DECISAOAOAOAO!=3 && DECISAOAOAOAO!=4 && DECISAOAOAOAO!=5 && DECISAOAOAOAO!=6)
    {
      cerr << "merd" << endl;
      exit(1);
    }

  if(DECISAOAOAOAO==1 || DECISAOAOAOAO==2 || DECISAOAOAOAO==3 || DECISAOAOAOAO==5)
    {
      cout << "Número complexo 2" << endl;
      cin >> c2;
      cout << endl;

      if(DECISAOAOAOAO==1)
	result=c1+c2;

      if(DECISAOAOAOAO==2)
	result=c1-c2;

      if(DECISAOAOAOAO==3)
	result=c1*c2;

      if(DECISAOAOAOAO==5)
	result=c1/c2; 
    }
  
  if(DECISAOAOAOAO==4 || DECISAOAOAOAO==6)
    {
      cout << "Escalar" << endl;;
      cin >> c;
      cout << endl;

      if(DECISAOAOAOAO==4)
	result=c*c1;

      if(DECISAOAOAOAO==6)
	result=c1/c;
    }

  cout << "Resultado" << endl;
  cout << result << endl;
     
  return 0;  

}

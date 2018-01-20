#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char const *argv[])
{
	int D, col, n;
	cout << "Dimensão: " << endl;
	cin >> D;
	cout << "Tamanho: " << endl;
	cin >> col;
	ofstream pot("potencial.txt");
	while(1)
	{
		cout << "Valor a colocar: " << endl;
		double v;
		cin >> v;
		int *coords = new int[D]();
		cout << "Indice a incrementar: " << endl;
		cin >> n;
		for(int i = 0; i<D; i++)
		{
			if(i != n)
			{
				cout << "Indice x" << i << ": " << endl;
				cin >> coords[i];
			}
		}

		coords[n] = 0;
		for (int i = 0; i<col; i++)
		{
			for(int j = 0; j<D; j++)
			{
				pot << coords[j] << "\t";
			}
			pot << v << endl;
			coords[n]++;
		}	

		int y;
		cout << "Mais? " << endl;
		cin >> y;

		if (y <= 0)
			break;
	}
	return 0;
}
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#define  PRES  double    



using namespace std;

double** createMatrix(int x)
{
	double** A = new double* [x];
	for (int i = 0; i < x; i++)
		A[i] = new double[x + 1];
	return A;
}
void deleteMatrix(double** X, int x)
{
	for (int i = 0; i < x; i++)
		delete X[i];
	delete[] X;
}
void initVector(double* X, int x, ifstream& in)
{
	for (int i = 0; i < x; i++)
	{
		double element; 
		in >> element;
		X[i] = element;
	}
}
void viewMatrix(double** X, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			cout << setw(15) << X[i][j];
		cout << endl;
	}
	cout << endl << endl;
}
void copyMatrix(double** X, double** copyX, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			copyX[i][j] = X[i][j];
	}
}
void zeroMatrix(double** X, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
		{
			X[i][j] = 0.0;
		}
	}
}
void viewAnswer(double* X, int x)
{
	for (int i = 0; i < x; i++)
		cout << "a[" << i << "] = " << X[i] << endl;
}
bool Gauss(double* An, double** X, int x)
{

	for (int k = 0; k < x; k++)
	{
		double max = fabs(X[k][k]);
		int remeber = k;		//запоминаем строку, чтобы не поменять саму себя
		for (int i = k + 1; i < x; i++)
		{
			if (max < fabs(X[i][k]))		//находим максимальный по модулю элемент в столбце
			{
				max = fabs(X[i][k]);
				remeber = i;
			}
		}

		if (fabs(max - 0) < 1e-6)
			return 0;

		if (k != remeber)				//меняем строки местами
		{
			double* temp = X[k];
			X[k] = X[remeber];
			X[remeber] = temp;
		}

		double lead = X[k][k];			//запоминаем ведущий элемент
		for (int r = k; r < x + 1; r++)
		{
			X[k][r] /= lead;
		}
		//начиная со следующей строки выполняем преобразование Гаусса
		for (int i = k + 1; i < x; i++)
		{
			double temp = X[i][k];
			for (int j = k; j < x + 1; j++)
			{
				X[i][j] -= X[k][j] * temp;
			}
		}
	}

	An[x - 1] = X[x - 1][x + 1 - 1];				//обратный ход
	for (int i = x - 2; i >= 0; i--)
	{
		An[i] = X[i][x + 1 - 1];
		for (int j = i + 1; j < x + 1 - 1; j++)
		{
			An[i] -= X[i][j] * An[j];
		}
	}
	return 1;
}

//	остаточная дисперсия S^2
double Dispersion(double* H, double* mu, double* An, int N, int m)
{
	double S = 0, temp;
	for (int i = 0; i < N; i++) {
		temp = mu[i];
		for (int j = 0; j < m + 1; j++)
			temp -= An[j] * pow(H[i], j);  //	по формуле 4.5
		S += temp * temp;
	}
	return S /= (N - m - 1);
}

void main()
{
	//variant 9
	//mu = a / H + b, a-?

	const int N = 6;	//	количество узлов
	int m = 2;	//	степень аппроксимирующего полинома
	double* H = new double[N];		//	напор жидкости
	double* mu = new double[N];		//	коэффициент истечения
	ifstream in_H("H.txt");
	ifstream in_mu("mu.txt");
	initVector(H, N, in_H);
	//a/H+b => a*x+b
	for (int i = 0; i < N; i++)
		H[i] = 1 / H[i];
	initVector(mu, N, in_mu);

	viewAnswer(H, N);
	viewAnswer(mu, N);

	double** A = createMatrix(m + 1);	
	zeroMatrix(A, m + 1);

	for (int i = 0; i < m + 1; i++)   //	составляем матрицу
	{
		for (int j = 0; j < m + 1; j++)
		{
			for (int k = 0; k < N; k++)
			{
				
				A[i][j] += pow(H[k], i + j);
			}
		}
	}

	for (int i = 0; i < m + 1; i++) //	столбец свободных членов
	{
		for (int j = 0; j < N; j++)
			A[i][m + 1] += mu[j] * pow(H[j], i);
	}

	viewMatrix(A, m + 1);

	double* An = new double[m + 1];
	double** copyA = createMatrix(m + 1);
	copyMatrix(A, copyA, m + 1);

	Gauss(An, copyA, m + 1);
	viewAnswer(An, m + 1);

	//	среднеквадратическое отклонение
	cout << "Sigma = " << sqrt(Dispersion(H, mu, An, N, m)) << endl;
//---------------------------------------------------------//
	int i, n;
	PRES w;

	ofstream foun("C:/Users/user/source/repos/chislaki/lab4/lab4/n.dat", ios_base::out | 
																		ios_base::trunc | 
																		ios_base::binary); 
	n = N;
	foun.write((char*)&n, sizeof n); 
	foun.close(); 

	PRES  x[N] = { 0.448, 0.432, 0.421, 0.417, 0.414, 0.412 };
	PRES  y1[N] = { 0.164, 0.328, 0.656, 0.984, 1.312, 1.640 };
	PRES  y2[N] = {};
	for (i = 0; i < N; i++) {
		y2[i] = An[i]*pow(x[i],i);
	}

	ofstream foutx("C:/Users/user/source/repos/chislaki/lab4/lab4/x.dat", ios_base::out | ios_base::trunc | ios_base::binary); 
	for (i = 0; i < N; i++) 
	{
		w = x[i];
		foutx.write((char*)&w, sizeof w);
	}

	for (i = 0; i < N; i++) 
	{
		w = y1[i];
		foutx.write((char*)&w, sizeof w);
	}

	for (i = 0; i < N; i++) 
	{
		w = y2[i];
		foutx.write((char*)&w, sizeof w);
	}

	foutx.close(); 

	delete An;
	delete H;
	delete mu;
	deleteMatrix(A, m + 1);
	deleteMatrix(copyA, m + 1);
}
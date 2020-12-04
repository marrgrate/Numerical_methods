#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

//��������� ������ � ����������� ���������
double** createMatrix(int x, int y)
{
	double** A = new double* [x];
	for (int i = 0; i < x; i++)
		A[i] = new double[y];
	return A;
}
//������������ ������
void deleteMatrix(double** X, int x)
{
	for (int i = 0; i < x; i++)
		delete X[i];
	delete[] X;
}
//���������� ������� ��������
void initMatrix(double** X, int rows, int columns, ifstream& in)
{
	for (int i = 0; i < rows; i++)
	{
		//cout << i + 1 << " equation\n";
		for (int j = 0; j < columns; j++)
		{
			double element; in >> element;
			X[i][j] = element;
		}
	}
}
//������ ������� ��������� ������ ����������� ������� �� ������� X �������
void reinitMatrix(double** A, double* X, int x, int y)
{
	for (int i = 0; i < x; i++)
	{
		A[i][y - 1] = 0;
		for (int j = 0; j < y - 1; j++)
		{
			A[i][y - 1] += A[i][j] * X[j];
		}
	}
}
//����������� � �������
void viewMatrix(double** X, int x, int y)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
			cout << setw(12) << X[i][j];
		cout << endl;
	}
	cout << endl << endl;
}
//����������� � ������� �������
void viewAnswer(double* X, int x)
{
	for (int i = 0; i < x; i++)
		cout << setw(12) << X[i];
	cout << endl;
}
//����� �������������
void findMaxElementOfVector(double* V, int size, double& max)
{
	max = fabs(V[0]);
	for (int i = 1; i < size; i++)
		if (max < fabs(V[i]))
			max = fabs(V[i]);
}
//����������� ��� ���������� ���������
void copyMatrix(double** X, double** copyX, int x, int y)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
			copyX[i][j] = X[i][j];
	}
}
//����� ������ � �������� ���
int Gauss(double* An, double** X, int row, int column)
{

	for (int k = 0; k < row; k++)
	{
		double max = abs(X[k][k]);
		int remember = k;		//���������� ������, ����� �� �������� ���� ����
		for (int i = k + 1; i < row; i++)
		{
			if (max < fabs(X[i][k]))		//������� ������������ �� ������ ������� � �������
			{
				max = fabs(X[i][k]);
				remember = i;
			}
		}

		if (fabs(max - 0) < 1e-6)
		{
			return 0;
		}

		if (k != remember)				//������ ������ �������
		{
			swap(X[k], X[remember]);
		}

		//viewMatrix(X, row, column);

		double lead = X[k][k];			//���������� ������� �������
		for (int r = k; r < column; r++)
		{
			X[k][r] /= lead;
		}
		//������� �� ��������� ������ ��������� �������������� ������
		for (int i = k + 1; i < row; i++)
		{
			double temp = X[i][k];
			for (int j = k; j < column; j++)
			{
				X[i][j] -= X[k][j] * temp;
			}
		}
		//viewMatrix(A, n, m);
	}

	An[row - 1] = X[row - 1][column - 1];				//�������� ���
	for (int i = row - 2; i >= 0; i--)
	{
		An[i] = X[i][column - 1];
		for (int j = i + 1; j < column - 1; j++)
		{
			An[i] -= X[i][j] * An[j];
		}
	}
	return 1;
}
//������ �������
void calculateNev(double* Nev, double* An, double** X, int row, int column)
{
	for (int i = 0; i < row; i++)
	{
		double sum = 0;
		for (int j = 0; j < column - 1; j++)
		{
			sum += X[i][j] * An[j];
		}
		Nev[i] = sum - X[i][column - 1];
	}
	double norma = 0;
	findMaxElementOfVector(Nev, row, norma);

	cout << "\n Norma = " << norma << endl << endl;
	cout << " Nevyazka\t";  viewAnswer(Nev, row);
}
//���������� �����������
void calculateDelta(double* XX, double* X, double* delta, int row)
{
	double* DIF = new double[row];
	for (int i = 0; i < row; i++)
		DIF[i] = XX[i] - X[i];
	cout << "\nDifference X~~ - X~\t"; viewAnswer(DIF, row);

	double x = 0, dif = 0;
	findMaxElementOfVector(DIF, row, dif); 
	cout << "\nmax |X~~ - X~|\t" << dif;
	findMaxElementOfVector(X, row, x);	
	cout << "\nmax |X~|\t" << x;
	*delta = dif / x;
}

void main()
{
	int row = 3; //cout << "\n Enter the number of equations: "; cin >> row;
	int column = 4; //cout << "\n Enter the number of roots of the equation: "; cin >> column; column++;

	double** A = createMatrix(row, column);
	
	ifstream in("mx.txt");
	initMatrix(A, row, column, in);
	in.close();

	double** B = createMatrix(row, column);
	copyMatrix(A, B, row, column);
	cout << endl << endl;
	viewMatrix(B, row, column);

	double* X = new double[row];

	if (Gauss(X, B, row, column) == 1)
	{
		viewMatrix(B, row, column);
		cout << endl << " Answer X~";
		viewAnswer(X, row);

		double* nevyazka = new double[row];
		calculateNev(nevyazka, X, A, row, column);
		delete nevyazka;
		
		reinitMatrix(A, X, row, column);
		cout << endl << endl;
		viewMatrix(A, row, column);

		double* XX = new double[row];	

		//���������� x~~
		if (Gauss(XX, A, row, column) == 1)
		{
			viewMatrix(A, row, column);
			cout << endl << "Answer X~~";  viewAnswer(XX, row);
			double delta;
			calculateDelta(XX, X, &delta, row);
			cout << "\nDelta = " << delta << endl;
		}
	}
	else cout << "\n The system has no solutions!\n\n";
	delete X;
	deleteMatrix(B, row);
	deleteMatrix(A, row);
}
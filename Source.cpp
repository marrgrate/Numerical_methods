#include <iostream>
#include <iomanip>


using namespace std;
const double de = 1e-9;
const int itr_max = 100;


double f1(double x1, double x2) {
	return 1.5 * pow(x1, 3) - pow(x2, 3) - 1;
}
double f2(double x1, double x2) {
	return x1 * pow(x2, 3) - x2 - 4;
}

double** createMatrix(int x)
{
	double** A = new double* [x];
	for (int i = 0; i < x; i++)
		A[i] = new double[x + 1];
	return A;
}
void deleteMatrix(double** X, int size)
{
	for (int i = 0; i < size; i++)
		delete X[i];
	delete[] X;
}
void copyMatrix(double** X, double** copyX, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			copyX[i][j] = X[i][j];
	}
}
void viewMatrix(double** X, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			cout << setw(12) << X[i][j]; 
		cout << endl;
	}
	cout << endl << endl;
}
void viewAnswer(double* X, int x)
{
	for (int i = 0; i < x; i++)
		cout << setw(12) << X[i];
	cout << endl;
}
bool Gauss(double* An, double** X, int x)
{
	for (int k = 0; k < x; k++)
	{
		double max = fabs(X[k][k]);
		int remeber = k;		//запоминаем строку, чтобы не помен€ть саму себ€
		for (int i = k + 1; i < x; i++)
		{
			if (max < fabs(X[i][k]))		//находим максимальный по модулю элемент в столбце
			{
				max = fabs(X[i][k]);
				remeber = i;
			}
		}

		if (fabs(max - 0) < 1e-6)
		{
			return 0;
		}

		if (k != remeber)				//мен€ем строки местами
		{
			double* temp = X[k];
			X[k] = X[remeber];
			X[remeber] = temp;
		}

		//viewMatrix(X, x);

		double lead = X[k][k];			//запоминаем ведущий элемент
		for (int r = k; r < x + 1; r++)
		{
			X[k][r] /= lead;
		}
		//начина€ со следующей строки приводим исходную матрицу к диагональному виду
		for (int i = k + 1; i < x; i++)
		{
			double temp = X[i][k];
			for (int j = k; j < x + 1; j++)
			{
				X[i][j] -= X[k][j] * temp;
			}
		}
		//viewMatrix(X, x);
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

typedef double(*pf)(double, double);

//подсчет производной конечно-разностным методом
double Differential(pf f, double x1, double x2, int x){
	if (x == 1)											//по первой неизвестной
		return (f(x1 + de, x2) - f(x1, x2)) / de;
	else                                                //по второй неизвестной
		return (f(x1, x2 + de) - f(x1, x2)) / de;
}

double Derivative1(double x1, double x2, int x){
	if (x == 1)											//по первой неизвестной
		return 4.5 * pow(x1, 2);
	else                                                //по второй неизвестной
		return -4 * x2;
}

double Derivative2(double x1, double x2, int x) {
	if (x == 1)											//по первой неизвестной
		return pow(x2, 3);
	else                                                //по второй неизвестной
		return 3 * x1 * pow(x2, 2) - 1;
}

//метод Ќьютона
double Newton(double* function, double* approximation, int n, int c) { //c задает метод подсчета производных
	double** F = createMatrix(n);
	double* dX = new double[n];	//delta x
	double e = 1e-9, delta1 = 0, delta2 = 0;

	cout << endl << setw(15) << "k" << setw(20) << "delta1      " << setw(15) << "delta2      " << endl << endl;
	int itr = 1;

	do
	{
		if (c == 1)
		{
			//якобиан, аналитический способ расчета производных
			F[0][0] = Derivative1(approximation[0], approximation[1], 1);
			F[0][1] = Derivative1(approximation[0], approximation[1], 2);
			F[0][2] = -f1(approximation[0], approximation[1]);					
			F[1][0] = Derivative2(approximation[0], approximation[1], 1);
			F[1][1] = Derivative2(approximation[0], approximation[1], 2);
			F[1][2] = -f2(approximation[0], approximation[1]);					
		}
		else
		{
			//заполнение матрицы якоби в €чейки F[0][0], F[0][1], F[1][0], F[1][1] конечно разностным методом
			F[0][0] = Differential(f1, approximation[0], approximation[1], 1);
			F[0][1] = Differential(f1, approximation[0], approximation[1], 2);
			F[0][2] = -f1(approximation[0], approximation[1]);					
			F[1][0] = Differential(f2, approximation[0], approximation[1], 1);
			F[1][1] = Differential(f2, approximation[0], approximation[1], 2);
			F[1][2] = -f2(approximation[0], approximation[1]);					
		}
		
		double** copyF = createMatrix(n);					// копи€ дл€ сохранени€ оригинала
		copyMatrix(F, copyF, n);

		if (!Gauss(dX, copyF, n))				
		{
			deleteMatrix(copyF, n);		
			deleteMatrix(F, n);
			delete dX;
			return 0;
		}

		approximation[0] += dX[0];				//уточнение решени€
		approximation[1] += dX[1];

		if (fabs(f1(approximation[0], approximation[1])) > fabs(f2(approximation[0], approximation[1]))) 
			//подсчет первой погрешности
			delta1 = fabs(f1(approximation[0], approximation[1]));
		else
			delta1 = fabs(f2(approximation[0], approximation[1]));

		for (int i = 0; i < n; i++)										
			//подсчет второй погрешности
		{
			if (fabs(dX[i]) < 1)
				delta2 = fabs(dX[i]);
			else if (fabs(dX[i]) >= 1)
				delta2 = fabs(dX[i] / approximation[i]);
		}

		cout << setw(15) << itr << setw(15) << delta1 << setw(15) << delta2 << endl;
		itr++;
	} while ((delta1 > e || delta2 > e) && (itr < itr_max));

	function[n - 2] = approximation[0];
	function[n - 1] = approximation[1];

	deleteMatrix(F, n);	
	delete dX;
	return *function;	//возвращаем решение —ЌЋј”
}

void main()
{
	//Var 9
	const int n = 2;
	double* x0 = new double[n];
	double* F = new double[n];
	x0[0] = 1;
	x0[1] = 1;

	*F = Newton(F, x0, n, 1);
	if (*F != 0)
	{
		cout << "\n\n Answer: ";
		viewAnswer(F, n);
	}
	else
		cout << "\n\n There is no solution!\n\n";

	x0[0] = 1;
	x0[1] = 1;
	*F = Newton(F, x0, n, 2);
	if (*F != 0)
	{
		cout << "\n\n Answer: ";
		viewAnswer(F, n);
	}
	else
		cout << "\n\n There is no solution!\n\n";
}
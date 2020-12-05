#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double Function1(double x)
{
	return (1 + x + x * x)/ sqrt(x * x * x - 1);
}

double Function2(double x, double y)
{
	return 1 / pow(x + y, 2);
}

double TrapezoidalRule(double a, double b, int& iter)
{
	double result = 0, prevResult;
	double h = b - a;

	do {
		prevResult = result;
		result = 0;
		for (int i = 1; i <= (b - a) / h - 1; ++i) 
		{
			result += 2 * Function1(a + h * i);
		}
		result += Function1(a) + Function1(b);
		result *= h / 2;
		h /= 2;
		iter++;
	} while (abs(result - prevResult) > 3 * 1e-4);
	return result;
}

double SimpsonsRule(double a, double b, int& count)
{
	double result = 0, result_prev = 0;
	double h = (b - a) / 2;
	do
	{
		result_prev = result;
		result = 0;
		for (int i = 1; i <= (b - a) / h; i += 2)
		{
			result += 4 * Function1(a + h * i);

		}
		for (int i = 2; i < (b - a) / h - 1; i += 2) {
			result += 2 * Function1(a + h * i);
		}
		result += Function1(a) + Function1(b);
		result *= h / 3;
		h /= 2;
		count++;
	} while (abs(result - result_prev) > 15 * 1e-4);
	return result;
}

double CubatureSimpsonsRule(double a, double b, double c, double d)
{
	int m = 2, n = 2 * m;
	double hx = (b - a) / (2 * n), hy = (d - c) / n;
	double result = 0.0;

	double xi = a;
	double yi = c;

	double* Xi = new double[2 * n + 1];
	Xi[0] = xi;

	for (int i = 1; i <= 2 * n; i++)
		Xi[i] = Xi[i - 1] + hx;;

	double* Yi = new double[2 * m + 1];
	Yi[0] = yi;

	for (int j = 1; j <= 2 * m; j++)
		Yi[j] = Yi[j - 1] + hy;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result += Function2(Xi[2 * i], Yi[2 * j]);
			result += 4 * Function2(Xi[2 * i + 1], Yi[2 * j]);
			result += Function2(Xi[2 * i + 2], Yi[2 * j]);
			result += 4 * Function2(Xi[2 * i], Yi[2 * j + 1]);
			result += 16 * Function2(Xi[2 * i + 1], Yi[2 * j + 1]);
			result += 4 * Function2(Xi[2 * i + 2], Yi[2 * j + 1]);
			result += Function2(Xi[2 * i], Yi[2 * j + 2]);
			result += 4 * Function2(Xi[2 * i + 1], Yi[2 * j + 2]);
			result += Function2(Xi[2 * i + 2], Yi[2 * j + 2]);
		}
	}
	result *= (hx * hy / 9);
	return result;
}


void main()
{
	//varint 9
	double a = 1.1, b = 2.631;

	cout << "\tTRAPEZOIDAL RULE: \n"
		<< "Intervals: a = 1.0 , b = 2.631\n";
	int iter = 0;	//количество итераций
	double result = TrapezoidalRule(a, b, iter);
	cout << " Result = " << result << " with " << iter << " iterations \n";


	cout << " \n\tSIMPSON'S RULE:\n"
		<< "Intervals: a = 1.0 , b = 2.631\n";
	result = 0;
	iter = 0;
	result = SimpsonsRule(a, b, iter);
	cout << " Result =  " << result << " with " << iter << " iterations \n";

	cout << " \n\tCUBATURE SIMPSON'S RULE:\n"
			<< "Intervals: a = 3.0 , b = 4.0 | c = 1.0 d = 2.0 \n";
	a = 3.0, b = 4.0, result = 0;
	double c = 1.0, d = 2.0;
	iter = 0;
	result = CubatureSimpsonsRule(a, b, c, d);
	cout << " Result = " << result << " \n";
}
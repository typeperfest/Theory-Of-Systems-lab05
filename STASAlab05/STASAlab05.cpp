#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <Windows.h>

using namespace std;

const int c = -4;
const int d = 2;
const int a = -1;
const int b = 2;
const int N = 12;
const int A = 5;
double c_min = -4.5, c_max = -3.5, d_min = 1.5, d_max = 2.5;
vector<double> X, Y, MY;
pair<double, double> weight;

void Param() {
	//	srand(static_cast<unsigned int>(time()));
	srand(3);
	for (double i = a; i < b; i += (b - a) * 1.0 / (N - 1)) {
		X.push_back(i);
		Y.push_back(i * c + d);
		MY.push_back(i * c + d + A * ((rand() % 10) / 10.0 - 0.5));
	}
}

double SKO(double c_, double d_, int f) {
	double sum = 0;
	if (f == 0) {
		for (size_t i = 0; i < N; ++i) {
			sum += pow(X[i] * c_ + d_ - MY[i], 2);
		}
	}
	else {
		for (size_t i = 0; i < N; ++i) {
			sum += pow(X[i] * c_ + d_ - Y[i], 2);
		}
	}
	return sum;
}

pair<double, double>MNK() {
	double sum_x = 0, sum_t = 0, sum_x_t = 0;
	for (size_t i = 0; i < N; ++i) {
		sum_x += X[i];
		sum_t += MY[i];
		sum_x_t += X[i] * MY[i];
	}
	double c_ = (N * sum_x_t - sum_x * sum_t) / (N * (sum_x * sum_x) - sum_x * sum_x);
	double d_ = (sum_t - c_ * sum_x) / N;
	return make_pair(c_, d_);
}

int F(int n) {
	int f, f1(1), f2(1), m(0);
	while (m < n - 1)
	{
		f = f1 + f2;
		f1 = f2;
		f2 = f;
		++m;
	}
	return f1;
}

double search_d(double c_, int f) {
	const double t = 1.618034, epsilon = 0.01;
	double myepsilon, ak = d_min, bk = d_max, x1, x2, y1, y2;
	int G = 100;
	x1 = ak + (double)F(G - 2) / F(G) * (bk - ak);
	x2 = ak + (double)F(G - 1) / F(G) * (bk - ak);
	if (f == 0) {
		y1 = SKO(c_, x1, 0);
		y2 = SKO(c_, x2, 0);
	}
	if (f == 1) {
		y1 = SKO(c_, x1, 1);
		y2 = SKO(c_, x2, 1);

	}
	for (; G <= 1; --G) {
		if (y1 > y2) {
			ak = x1;
			x1 = x2;
			x2 = bk - (x2 - ak);
			y1 = y2;
			y2 = SKO(c_, x2, 0);
		}
		else {
			bk = x2;
			x2 = x1;
			x1 = ak + (bk - x2);
			y2 = y1;
			y1 = SKO(c_, x1, 0);
		}
	}
	return (x1 + x2) / 2;
}

pair<double, double> search_c(int f) {
	const double epsilon = 0.01;
	double myepsilon, ak = c_min, bk = c_max, x, delta, min_c = 2;
	for (size_t i = 1; i <= N; ++i) {
		delta = ((b - a) / (i + 1.));
		if (f == 0) {
			for (size_t j = 1; j <= i; ++j) {
				x = ak + (j * delta);
				if (SKO(x, search_d(x, 0), 0) <= SKO(min_c, search_d(min_c, 0), 0)) {
					min_c = x;
				}
			}
		}
		else {
			for (size_t j = 1; j <= i; ++j) {
				x = ak + (j * delta);
				if (SKO(x, search_d(x, 1), 1) <= SKO(min_c, search_d(min_c, 1), 1)) {
					min_c = x;
				}
			}
		}
	}
	if (f == 0) {
		return make_pair(min_c, search_d(min_c, 0));
	}
	if (f == 1) {
		return make_pair(min_c, search_d(min_c, 1));
	}
}
void print(const vector<double> & v) {
	for (size_t i = 0; i < v.size(); ++i) {
		cout << setprecision(3) << i + 1 << ") " << '(' << X[i] << ',' << v[i] << ')' << endl;
	}
}

void print(const double& c_, const double& d_) {
	for (size_t i = 0; i < N; ++i) {
		cout << setprecision(3) << i + 1 << ") " << '(' << X[i] << ',' << X[i] * c_ + d_ << ')' << endl;
	}
}
int main()
{
	setlocale(LC_ALL, "RUSS");
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	Param();
	cout << "Множество значений y(x) с коэффициентами с=" << c << ", " << "d=" << d << ": " << endl;
	print(Y);
	cout << "Множество значений зашумленной функции y(x) с коэффициентами с=" << c << ", " << "d=" << d << ", заданных равномерно со случайными ошибками" << ':' << endl;
	print(MY);
	auto w = MNK();
	cout << "Множество значений функции f(x) с коэффициентами с*=" << w.first << ", " << "d*=" << w.second << ", найдеными по МНК" << ':' << endl;
	print(w.first, w.second);
	auto w1 = search_c(1);
	cout << "Множество значений функции g(x) с набором синаптических весов " << "c*=" << w1.first << ", " << "d*=" << w1.second << ",найденных для функции без шума :" << endl;
	print(w1.first, w1.second);
	auto w0 = search_c(0);
	cout << "Множество значений функции h(x) с набором синаптических весов " << "c*=" << w0.first << ", " << "d*=" << w0.second << ", найденных для функции c шумом :" << endl;
	print(w0.first, w0.second);
	system("pause");
	return 0;
}

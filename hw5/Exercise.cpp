#include "Exercise.h"
#include"Function.h"
#include<time.h>
#include<iostream>
#include <cmath>

void exercise1() {
	int n = 20;
	int N = (n - 1) * (n - 1);
	double h = static_cast<double>(1) / n;
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	clock_t start = 0, end = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				A[i][j] = 1 + h * h / 4;
			}
			else if ((j == i + 1 && j % (n - 1) != 0 )|| (j == i - 1 && i % (n - 1) != 0)) {
				A[i][j] = static_cast<double>(-1) / 4;
			}
			else if (j == i + n - 1 || j == i - n + 1) {
				A[i][j] = static_cast<double>(-1) / 4;
			}
			else {
				A[i][j] = 0;
			}
		}
	}
	//matrix_print(A);
	
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			b[i * (n - 1) + j] = h * h * sin((i + 1) * (j + 1) * h * h) / 4;
			if (i == 0) {
				b[i * (n - 1) + j] += h * h * (j + 1) * (j + 1) / 4;
			}
			if (j == 0) {
				b[i * (n - 1) + j] += h * h * (i + 1) * (i + 1) / 4;
			}
			if (i == n - 2) {
				b[i * (n - 1) + j] += static_cast<double>(1) / 4+ h * h * (j + 1) * (j + 1) / 4;
			}
			if (j == n - 2) {
				b[i * (n - 1) + j] += static_cast<double>(1) / 4+ h * h * (i + 1) * (i + 1) / 4;
			}
		}
	}
	start = clock();
	int k = conjugate_gradient(A, b);
	end = clock();
	cout << "共轭梯度法解向量为" << endl;
	solution_print(b);
	
	cout << "迭代次数为" << k << endl;
	cout << "运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	/*
	start = clock();
	k = sor(A, b,1);
	end = clock();
	cout << "SOR迭代法解向量为" << endl;
	solution_print(b);
	
	cout << "迭代次数为" << k << endl;
	cout << "运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;*/
}

void exercise2(int N) {
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	clock_t start=0, end=0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			A[i][j] = 1 / (static_cast<double>(i + j + 1));
		}
	}
	for (int i = 0; i < N; i++) {
		b[i] = 0;
		for (int j = 0; j < N; j++)
		{
			b[i] = b[i] + 1 / (static_cast<double>(i + j + 1));
		}
		b[i] = b[i] / 3;
	}
	start = clock();
	int k=conjugate_gradient(A, b);
	end = clock();
	cout << N << "阶Hilbert矩阵的共轭梯度法解向量为" << endl;
	solution_print(b);
	cout<<"迭代次数为" << k << endl;
	cout << "运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl<<endl;
}

void exercise3() {
	int N = 5;
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	clock_t start = 0, end = 0;
	A = { {10,1,2,3,4},
		{1,9,-1,2,-3},
		{2,-1,7,3,-5},
		{3,2,3,12,-1},
		{4,-3,-5,-1,15} };
	b = { 12,-27,14,-17,12 };
	//共轭梯度法
	start = clock();
	int k = conjugate_gradient(A, b);
	end = clock();
	cout << "共轭梯度法解向量为" << endl;
	solution_print(b);
	cout << "迭代次数为" << k << endl;
	cout << "运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;

	//Jacobi迭代
	
	b = { 12,-27,14,-17,12 };
	start = clock();
	k = jacobi(A, b);
	end = clock();
	cout << "Jacobi迭代法解向量为" << endl;
	solution_print(b);
	cout << "迭代次数为" << k << endl;
	cout << "运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;

	//G-S迭代
	
	b = { 12,-27,14,-17,12 };
	start = clock();
	k = g_s(A, b);
	end = clock();
	cout << "G-S迭代法解向量为" << endl;
	solution_print(b);
	cout << "迭代次数为" << k << endl;
	cout << "运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;

	//SOR迭代

	b = { 12,-27,14,-17,12 };
	start = clock();
	k = sor(A, b, 1.3);
	end = clock();
	cout << "SOR迭代法解向量为" << endl;
	solution_print(b);
	cout << "迭代次数为" << k << endl;
	cout << "运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;

}


void test() {/*
	double w = h * h / 4;
	for (int i = 0; i < N; i++) {
		if (i == 0) {
			b[0] = w * sin(h * h) + 0.5 * 4 * w;
			b[n - 1 - 1] = w * sin((n - 1) * h * h) + 0.25 * (1 + h * h + N * h * h);
			for (int j = 1; j < n - 1 - 1; j++) {
				b[j] = w * sin((j + 1) * h * h) + 0.25 * (j + 1) * h * (j + 1) * h;
			}
			i = n - 1 - 1;
		}
		else if (i == (n - 1) * (n - 1 - 1)) {
			b[((n - 1) * (n - 1) - 1)] = w * sin(h * (n - 1) * h) + 0.25 * (1 + h * h + N * h * h);
			b[(n - 1) * (n - 2)] = w * sin(N * h * h) + 0.25 * (2 + 2 * N * h * h);
			for (int j = 1; j < n - 2; j++) {
				b[(n - 1) * (n - 2) + j] = w * sin((j + 1) * (n - 1) * h * h) + 0.25 * (1 + (j + 1) * (j + 1) * h * h);
			}
			i = N;
		}
		else {
			int P = i % (n - 1) + 1, Q = i / (n - 1) + 1;
			if (P == 1) {
				b[i] = w * sin(P * Q * h * h) + 0.25 * Q * Q * h * h;
			}
			else if (P == n - 1) {
				b[i] = w * sin(P * Q * h * h) + 0.25 * (1 + Q * Q * h * h);
			}
			else {
				b[i] = w * sin(P * Q * h * h);
			}
		}
	}
*/
}
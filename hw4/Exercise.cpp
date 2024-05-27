#include "Exercise.h"
#include"Function.h"
#include<time.h>
#include<iostream>
#include <cmath>

int exercise1(double epsilon, double  omega) {
	double a = 0.5, n = 100.0, h = 1 / n;
	int k = 0;//iteration times
	clock_t start = 0, end = 0;
	vector<double> z(n-1), y(n-1), x(n-1), b(n-1);
	vector<vector<double>> A(n - 1, vector<double>(n - 1));
	for (int i = 0; i < n-1; i++) {
		x[i] = (i+1) * h;
		z[i] = a * x[i] + (1 - a) * (1 - exp(-x[i] / epsilon)) / 
			(1 - exp(-1 / epsilon));
	}
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			if (i == j) {
				A[i][j] = -2 * epsilon - h;
			}
			else if (j == i + 1) {
				A[i][j] = epsilon + h;
			}
			else if (j == i - 1) {
				A[i][j] = epsilon;
			}
			else {
				A[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < n - 1; i++) {
		b[i] = a * h * h;
	}
	b[n - 2] = a * h * h - epsilon - h;
	cout << "n="<<n<<",epsilon="<<epsilon<<",omega="<<omega<<"时" << endl<<endl;
	//solution_print(z);
	//Jacobi
	start = clock();
	k = jacobi(A, b);
	end = clock();
	for (int i = 0; i < n - 1; i++) {
		b[i] = b[i]-z[i];
	}
	cout << "the time of Jacobi " << double(end - start) / CLOCKS_PER_SEC << "s" ;
	cout << ", iteration times " << k << ",error "<<vector_infinity_norm(b)<<endl;

	//G-S
	
	for (int i = 0; i < n - 1; i++) {
		b[i] = a * h * h;
	}
	b[n - 2] = a * h * h - epsilon - h;
	start = clock();
	k = g_s(A, b);
	end = clock();
	for (int i = 0; i < n - 1; i++) {
		b[i] = b[i] - z[i];
	}
	cout << "G-S iteration running time " << double(end - start) / CLOCKS_PER_SEC << "s";
	cout << "iteration times " << k << ",error " << vector_infinity_norm(b) << endl;

	//SOR
	
	for (int i = 0; i < n - 1; i++) {
		b[i] = a * h * h;
	}
	b[n - 2] = a * h * h - epsilon - h;
	start = clock();
	k=sor(A,b,omega);
	end = clock();
	for (int i = 0; i < n - 1; i++) {
		b[i] = b[i] - z[i];
	}
	cout << "SOR iteration running time " << double(end - start) / CLOCKS_PER_SEC << "s";
	cout << " iteration times " << k << ",error " << vector_infinity_norm(b) << endl<<endl;
	return k;
}

int exercise2(int N, double omega) {
	cout << "N=" << N << ",omega=" << omega << "时，迭代结果为" << endl<<endl;
	vector<vector<double>> U(N, vector<double>(N)), F = U, G = U, W = U;
	                                 
	vector<double> b(N * N);         //表示x_k+1-x_k 和误差
	double h = static_cast<double>(1) / N;
	clock_t start = 0, end = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == 0 || j == 0) {
				U[i][j] = 1;
			}
			else {
				U[i][j] = 0;
			}
			b[i + j] = 1;
			F[i][j] = i * h + j * h;
			G[i][j] = exp(i * h * j * h);
		}
	}
	W = U;
	//Jacobi迭代
	int k = 0;
	start = clock();
	while (vector_2_norm(b) > 1e-7) {
		for (int i = 1; i < N-1; i++) {
			for (int j = 1; j < N-1; j++) {
				W[i][j] = (U[i - 1][j] + U[i + 1][j] + U[i][j + 1] + U[i][j - 1] +
					h * h * F[i][j]) / (4 + h * h * G[i][j]);
			}
		}
		for (int i = 1; i < N - 1; i++) {
			W[i][N - 1] = (U[i - 1][N - 1] + U[i + 1][N - 1] + 1 + U[i][N - 2] +
				h * h * F[i][N - 1]) / (4 + h * h * G[i][N - 1]);
		}
		for (int j = 1; j < N - 1; j++) {
			W[N-1][j] = (U[N-2][j] + 1 + U[N-1][j + 1] + U[N-1][j - 1] +
				h * h * F[N-1][j]) / (4 + h * h * G[N-1][j]);
		}
		W[N - 1][N - 1] = (U[N - 1][N - 2] + 1 + U[N - 2][N - 1] + 1 + 
			h * h * F[N - 1][N - 1])/ (4 + h * h * G[N - 1][N - 1]);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				b[i + j] = W[i][j] - U[i][j];
			}
		}
		U = W;
		k++;
	}
	end = clock();
	cout << "the time of Jacobi " << double(end - start) / CLOCKS_PER_SEC << "s";
	cout << "，iteration times为" << k << endl <<
		"解的最小分量" << matrix_min(U) << endl ;

	//G-S迭代
	//W表示U_k-1,即上一次的迭代结果
	k = 0;
	start = clock();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == 0 || j == 0) {
				U[i][j] = 1;
			}
			else {
				U[i][j] = 0;
			}
			b[i + j] = 1;
		}
	}
	W = U;
	while (vector_2_norm(b) > 1e-7) {
		for (int i = 1; i < N - 1; i++) {
			for (int j = 1; j < N - 1; j++) {
				U[i][j] = (U[i - 1][j] + U[i + 1][j] + U[i][j + 1] + U[i][j - 1] +
					h * h * F[i][j]) / (4 + h * h * G[i][j]);
			}
		}
		for (int i = 1; i < N - 1; i++) {
			U[i][N - 1] = (U[i - 1][N - 1] + U[i + 1][N - 1] + 1 + U[i][N - 2] +
				h * h * F[i][N - 1]) / (4 + h * h * G[i][N - 1]);
		}
		for (int j = 1; j < N - 1; j++) {
			U[N - 1][j] = (U[N - 2][j] + 1 + U[N - 1][j + 1] + U[N - 1][j - 1] +
				h * h * F[N - 1][j]) / (4 + h * h * G[N - 1][j]);
		}
		U[N - 1][N - 1] = (U[N - 1][N - 2] + 1 + U[N - 2][N - 1] + 1 +
			h * h * F[N - 1][N - 1]) / (4 + h * h * G[N - 1][N - 1]);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				b[i + j] = W[i][j] - U[i][j];
			}
		}
		W = U;
		k++;
	}
	end = clock();
	cout << "G-S iteration running time " << double(end - start) / CLOCKS_PER_SEC << "s";
	cout << "，iteration times为" << k << endl <<
		"解的最小分量" << matrix_min(U) << endl ;

	//SOR迭代
	//W表示U_k-1,即上一次的迭代结果
	k = 0;
	start = clock();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == 0 || j == 0) {
				U[i][j] = 1;
			}
			else {
				U[i][j] = 0;
			}
			b[i + j] = 1;
		}
	}
	W = U;
	while (vector_2_norm(b) > 1e-7) {
		for (int i = 1; i < N - 1; i++) {
			for (int j = 1; j < N - 1; j++) {
				U[i][j] = (U[i - 1][j] + U[i + 1][j] + U[i][j + 1] + U[i][j - 1] +
					h * h * F[i][j]) / (4 + h * h * G[i][j]);
				U[i][j] = W[i][j] + omega * (U[i][j] - W[i][j]);
			}
		}
		for (int i = 1; i < N - 1; i++) {
			U[i][N - 1] = (U[i - 1][N - 1] + U[i + 1][N - 1] + 1 + U[i][N - 2] +
				h * h * F[i][N - 1]) / (4 + h * h * G[i][N - 1]);
			U[i][N - 1] = W[i][N - 1] + omega * (U[i][N - 1] - W[i][N - 1]);
		}
		for (int j = 1; j < N - 1; j++) {
			U[N - 1][j] = (U[N - 2][j] + 1 + U[N - 1][j + 1] + U[N - 1][j - 1] +
				h * h * F[N - 1][j]) / (4 + h * h * G[N - 1][j]);
			U[N - 1][j] = W[N - 1][j] + omega * (U[N - 1][j] - W[N - 1][j]);
		}
		U[N - 1][N - 1] = (U[N - 1][N - 2] + 1 + U[N - 2][N - 1] + 1 +
			h * h * F[N - 1][N - 1]) / (4 + h * h * G[N - 1][N - 1]);
		U[N - 1][N - 1] = W[N - 1][N - 1] + omega * (U[N - 1][N - 1] - W[N - 1][N - 1]);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				b[i + j] = W[i][j] - U[i][j];
			}
		}
		W = U;
		k++;
	}
	end = clock();
	cout << "SOR iteration running time " << double(end - start) / CLOCKS_PER_SEC << "s";
	cout << "，iteration times为" << k << endl <<
		"解的最小分量" << matrix_min(U) << endl << endl;
	return k;
}

void test() {
	double omega = 0.1;
	int k = 0;
	k = exercise2(60, 0.1);
	for (double t = 1.896; t < 1.898; t += 0.0001) {
		if (k > exercise2(60, t)) {
			k = exercise2(60, t);
			omega = t;
		}
	}
	cout << omega;
}
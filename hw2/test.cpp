#include <iostream>
#include <vector>
#include<Windows.h>
#include"Function.h"
#include"Exercise.h"
#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<iostream>
using namespace std;

int main() {
	int N = 30;
	while (N >= 5) {
	int i = 0, j = 0;
	double k = 0.0, p = 0,p1=0;
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N), x(N), t(N);
	vector<int> u(N);
	

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++) {
				A[i][j] = 1 / (static_cast<double>(i + j + 1));
			}
		}

		for (int i = 0; i < N; i++) {
			b[i] = x[i] = 1;
		}
		
		matrix_vector_multiply(A, b);//Ax=b
		p1 = matrix_infinity_norm(A);
		p = vector_infinity_norm(b);
		gauss_elim_col_pivoting(A, u);
		p1 = p1 * matrix_inverse_infinity_norm(A);
		vector_pb(u, b);
		forward_subs1(A, b);
		back_subs(A, b);               //此时b是计算解

		for (int i = 0; i < N; i++) {
			t[i] = abs(b[i] - x[i]);
		}

		k = vector_infinity_norm(t) / vector_infinity_norm(x);

		cout << N << "阶该方阵的真实相对误差为" << k << "," ;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++) {
				A[i][j] = 1 / (static_cast<double>(i + j + 1));
			}
		}//重新赋值A，以计算Ax-Ax^
		
		matrix_vector_multiply(A, b);
		
		matrix_vector_multiply(A, x);
		
		for (int i = 0; i < N; i++) {
			t[i] = b[i] - x[i];
		}

		
		k = p1 * vector_infinity_norm(t)
			/ p;
		cout << "解的精度为" << k << endl;
		N--;
	}
	return 0;
}
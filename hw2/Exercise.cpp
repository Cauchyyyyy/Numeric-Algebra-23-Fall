#include "Exercise.h"
#include"Function.h"
#include<time.h>
#include<iostream>

void exercise_1(int N) {
	vector<vector<double>> A(N, vector<double>(N));
	vector<int> u(N);
	double k = 0;
	//输入Hilbert矩阵的值
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = 1 / (double(i + j + 1));
		}
	}
	k = matrix_infinity_norm(A);
	gauss_elim_col_pivoting(A, u);
	k = k* matrix_inverse_infinity_norm(A);
	cout << k << endl;
}

void exercise_2(int N) {
	int i = 0, j = 0;
	double k = 0.0,p=0,p1=0;
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N),x(N),t(N);
	vector<int> u(N);
	
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i == j || j == N - 1)
				A[i][j] = 1;
			else if (i > j)
				A[i][j] = -1;
			else
				A[i][j] = 0;
		}
	}
	srand((unsigned)time(NULL));  //x的随机赋值(来自网络)
	for (int i = 0; i < N; i++) {
		b[i] = x[i] = rand();
	}
	cout << "此时x的取值为" << endl;
	solution_print(x);

	
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
	
	k = vector_infinity_norm(t)/vector_infinity_norm(x);
	cout << N << "阶该方阵的真实相对误差为" << k<<",";
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i == j || j == N - 1)
				A[i][j] = 1;
			else if (i > j)
				A[i][j] = -1;
			else
				A[i][j] = 0;
		}
	}//重新赋值A，以计算Ax-Ax^
	
	matrix_vector_multiply(A, b);
	
	matrix_vector_multiply(A, x);
	
	for (int i = 0; i < N; i++) {
		t[i] = b[i] - x[i];
	}
	
	
	k=p1* vector_infinity_norm(t)
		/ p;
	cout <<"解的精度为" << k << endl;
}
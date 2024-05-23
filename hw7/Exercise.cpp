#include "Exercise.h"
#include"Function.h"
#include<time.h>
#include<iostream>
#include <cmath>

void exercise1(int n) {
	cout << "�������n=" << n << endl;
	vector<vector<double>> A(n, vector<double>(n));
	vector<vector<double>> Q = A;
	vector<double> ev(n);//����ֵ���У���������
	clock_t start = 0, end = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (abs(j - i) == 1) {
				A[i][j] = 1;
			}
			else {
				A[i][j] = 0;
			}
		}
		A[i][i] = 4;
	}
	start = clock();
	//jacobi_method2(A,Q); cout << "A_k=" << endl; matrix_print(A); cout << "Q_k=" << endl; matrix_print(Q);
	jacobi_method(A);
	end = clock();
	cout << "����ʱ��Ϊ" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << "�����������������ֵΪ" << endl;
	for (int i = 0; i < n; i++) {
		ev[i] = A[i][i];
	}sort_s_to_l(ev);
	solution_print(ev);
	
	cout << endl<<endl;
}

void exercise2() {
	int n = 100;
	vector<vector<double>> A(100, vector<double>(100));
	clock_t start = 0, end = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (abs(j - i) == 1) {
				A[i][j] = -1;
			}
			else { 
				A[i][j] = 0;
			}
		}
		A[i][i] = 2;
	}start = clock();
	double t = bisection(A, 1);
	end = clock();cout << "��С����ֵΪ" << t << endl;
	cout << "����ʱ��Ϊ" << double(end - start) / CLOCKS_PER_SEC << "s" << endl<<endl ;
	vector<double> u(n); 
	for (int i = 0; i < n; i++) {
		u[i] = 1;
	}cout << "��Ӧ����������";
	inverse_power_method(A, t, u,1); 
	//inverse_power_method2(A, t, u); 
	solution_print(u);

	cout << endl;
	start = clock();
	t = bisection(A, 100);
	end = clock();cout << "�������ֵΪ" << t << endl;
	cout << "����ʱ��Ϊ" << double(end - start) / CLOCKS_PER_SEC << "s" << endl<<endl;
	for (int i = 0; i < n; i++) {
		u[i] = 2;
	}cout << "��Ӧ����������";
	inverse_power_method(A, t, u,3); 
	//inverse_power_method2(A, t, u);
	solution_print(u);
}

void test() {
	
}
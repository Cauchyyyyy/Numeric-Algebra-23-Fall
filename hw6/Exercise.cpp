#include "Exercise.h"
#include"Function.h"
#include<time.h>
#include<iostream>
#include <cmath>

void exercise1() {
	vector<double> a,u;
	int n = 0;
	double x = 0;
	a = { 3,-5,1 };
	n = a.size();
	clock_t start = 0, end = 0;
	cout << "����ʽΪx^"<<n ;
	for (int i = n - 1; i >=0; i--) {
		if (a[i] != 0) {
			if (a[i] > 0) {
				cout << "+";
			}
			cout << a[i] << "x^" << i;
		}
	}cout << endl;
	start = clock();
	x = max_abs_root_monic_polynomial(a,a);
	end = clock();
	cout << "ģ����Ϊ"<<x << endl;
	cout << "����ʱ��Ϊ" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;
	a = {-1,-3,0};
	n = a.size();
	cout << "����ʽΪx^" << n;
	for (int i = n - 1; i >= 0; i--) {
		if (a[i] != 0) {
			if (a[i] > 0) {
				cout << "+";
			}
			cout << a[i] << "x^" << i;
		}
	}cout << endl;
	start = clock();
	x = max_abs_root_monic_polynomial(a, a);
	end = clock();
	cout << "ģ����Ϊ" << x << endl ;
	cout << "����ʱ��Ϊ" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;
	a = {-1000,790,-99902,79108.9,9802.08,10891.01,208.01,101};
	n = a.size();
	cout << "����ʽΪx^" << n;
	for (int i = n - 1; i >= 0; i--) {
		if (a[i] != 0) {
			if (a[i] > 0) {
				cout << "+";
			}
			cout << a[i] << "x^" << i;
		}
	}cout << endl;
	start = clock();
	x = max_abs_root_monic_polynomial(a, a);
	end = clock();
	cout << "ģ����Ϊ" << x << endl ;
	cout << "����ʱ��Ϊ" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;
}

void exercise2() {
	int n = 41;
	int k = 0;
	double x = 0.9;
	clock_t start = 0, end = 0;
	vector<double> a(41);
	a[3] = a[0] = 1;
	
	vector<vector<double>> A(n, vector<double>(n)),Q=A;//�����ѷ���
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == 0) {
				A[i][j] = -a[n - 1 - j];
			}
			else if (i == j + 1) {
				A[i][j] = 1;
			}
			else {
				A[i][j] = 0;
			}
		}
	}
	start = clock();
	k=implicit_qr_algorithm(A, Q); 
	end = clock();
	//cout << "��ʽQR�ֽ�õ���ʵSchur��׼�;���Ϊ" << endl;matrix_print(A); 
	cout << "�ö���ʽ�ĸ�Ϊ" << endl;
	find_eigenvalue_from_schur(A);
	cout<<"��������Ϊ"<<k<<"��,����ʱ��" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	vector<vector<double>> B(4, vector<double>(4)),P=B;
	
	while (x < 1.2) {
		B = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,x},{6.1,4.9,3.5,6.2} };
		start = clock();
		k = implicit_qr_algorithm(B, P);
		end = clock();
		cout <<"x="<<x << "ʱ����ʽQR�ֽ�õ���ʵSchur��׼�;���Ϊ" << endl;
		matrix_print(B);
		cout << "�þ��������ֵΪ" << endl;
		find_eigenvalue_from_schur(B);
		cout << "��������Ϊ" << k << "��,����ʱ��" << double(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;
		x += 0.1;
	}

}

void test() {
	int n = 4;
	
	vector<vector<double>> A(n, vector<double>(n)),Q=A;
	
	A = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,0.9},{6.1,4.9,3.5,6.2} };

	implicit_qr_algorithm(A, Q); matrix_print(A); find_eigenvalue_from_schur(A);
	//up_hessenberg_decomposition(A, Q);matrix_print(A);
	//matrix_multiply2(A, Q);matrix_full_transpose(Q); matrix_multiply(Q, A); matrix_print(A);
		
	//direct_qr_algorithm(A, Q, 30); matrix_print(A); find_eigenvalue_from_schur(A);
	//francis_displacement(A, Q);matrix_print(A);

	//matrix_multiply2(A, Q);matrix_full_transpose(Q); matrix_multiply(Q, A); matrix_print(A);
	
}
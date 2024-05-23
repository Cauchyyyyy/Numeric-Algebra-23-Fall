#include "Exercise.h"
#include"Function.h"
#include<time.h>
#include<iostream>

void exercise1_1() {
	int N = 84; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N), d(N);
	vector<int>u(N), v(N);
	double max;
	clock_t start, end;      //定义clock_t变量，用于计算算法所需时间

	//初始化A和b,QR分解法求解

	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 6;
		A[i + 1][i] = 8;
		A[i][i + 1] = 1;
		b[i] = 15;
	}
	A[N - 1][N - 1] = 6;
	b[0] = 7;
	b[N - 1] = 14;
	start = clock();           //运算开始时间（只计算算法用时，输出解等步骤不计入）
	qr_decomposition2(A, b, d);
	
	back_subs(A, b);
	end = clock();        //运算结束
	solution_print(b);
	cout << "求解矩阵1运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	//输出时间（单位：ｓ），这条语句来自网络，用于计算算法运行用时
	// 但在时间很短时显示为0s

	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout << "与精确解的最大误差为(一般采用x-x^的无穷范数表示)" << max << endl;
	cout << endl;

	for (int i = 0; i < N; i++)
	{
		b[i] = 15;
		for (int j = 0; j < N; j++) {
			if (j == i + 1)
				A[i][j] = 1;
			else if (j == i - 1)
				A[i][j] = 8;
			else if (j == i)
				A[i][j] = 6;
			else
				A[i][j] = 0;
		}
	}
	b[0] = 7;
	b[N - 1] = 14;
	
	start = clock();
	gauss_elim_col_pivoting(A, u);//与列主元gauss消去法的对比
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();
	
	cout << "列主元Gauss消去法运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout << "与精确解的最大误差为" << max << endl<<endl;
}

void exercise1_2() {
	int N = 100; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N), d(N), t(N);
	vector<int>u(N);
	double max;
	clock_t start, end;

	//初始化A和b
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			if (j == i)
				A[i][j] = 10;
			else if (j == i + 1 || j == i - 1)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}
	}

	srand((unsigned)time(NULL));  //b的随机赋值(来自网络)
	for (int i = 0; i < N; i++) {
		b[i] = t[i] = rand();
	}
	cout << "随机向量b的取值此时为" << endl;
	solution_print(b);

	cout << "相应的解如下" << endl << endl;

	start = clock();
	qr_decomposition(A, b, d);

	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "求解矩阵2运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			if (j == i)
				A[i][j] = 10;
			else if (j == i + 1 || j == i - 1)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}

	}
	matrix_vector_multiply(A, b);
	max = abs(b[0] - t[0]);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - t[i]))
			max = abs(t[i] - b[i]);
	}
	cout << "与精确解的最大误差为(这里采用的是Ax-Ab的无穷范数表示)" << max << endl;
	cout << endl;
	for (int i = 0; i < N; i++) {
		b[i] = t[i] ;
	}
	start = clock();
	modified_cholesky_decomp(A);
	forward_subs1(A, b);
	matrix_DLT(A);
	back_subs(A, b);
	end = clock();
	
	cout << "改进的平方根法运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			if (j == i)
				A[i][j] = 10;
			else if (j == i + 1 || j == i - 1)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}
	}
	matrix_vector_multiply(A, b);
	max = abs(b[0] - t[0]);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - t[i]))
			max = abs(t[i] - b[i]);
	}
	cout << "与精确解的最大误差为(这里采用的是Ax-Ab的无穷范数表示)" << max << endl;
	cout << endl;

	for (int i = 0; i < N; i++) {
		b[i] = t[i] ;
	}
	start = clock();
	gauss_elim_col_pivoting(A, u);
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();
	
	cout << "列主元Gauss消去法运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			if (j == i)
				A[i][j] = 10;
			else if (j == i + 1 || j == i - 1)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}
	}
	matrix_vector_multiply(A, b);
	max = abs(b[0] - t[0]);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - t[i]))
			max = abs(t[i] - b[i]);
	}
	cout << "与精确解的最大误差为(这里采用的是Ax-Ab的无穷范数表示)" << max << endl;
	cout << endl;
}

void exercise1_3(){
	int N = 10; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N), d(N);
	vector<int>u(N);
	double max;
	clock_t start, end;
	//初始化A和b
	
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
	}
	start = clock();
	qr_decomposition(A, b, d);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "求解Hilbert矩阵运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl ;

	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout << "与精确解的最大误差为" << max << endl;
	cout << endl;

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
	}
	start = clock();
	gauss_elim_col_pivoting(A, u);//与列主元gauss消去法的对比
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();

	cout << "列主元Gauss消去法运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout << "与精确解的最大误差为" << max << endl << endl;

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
	}
	start = clock();
	modified_cholesky_decomp(A);
	forward_subs1(A, b);
	matrix_DLT(A);
	back_subs(A, b);
	end = clock();
	cout << "改进的平方根法运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout << "与精确解的最大误差为" << max << endl << endl;
}

void exercise_2() {
	vector<double> t = { -1,-0.75,-0.5,0,0.25,0.5,0.75 };
	vector<double> y = { 1,0.8125,0.75,1,1.3125,1.75,2.3125 }, d = y;
	double r = 0;
	clock_t start = 0, end = 0;
	int m = t.size();
	vector<vector<double>> A(m, vector<double>(3));
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < 3; j++) {
			A[i][j] = pow(t[i], j);
		}
	}
	start = clock();
	qr_decomposition(A, y, d);
	back_subs(A, y);
	end = clock();
	cout << "a,b,c的值分别为" << y[2] <<" " << y[1] <<" "<< y[0] << endl;

	cout << "求解最小二乘运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	for (int i = 3; i < m; i++) {
		r = r + y[i] * y[i];
	}
	r = sqrt(r);
	cout << "残向量的2范数为" << r << endl << endl;

}

void exercise_3() {
	double r = 0;
	clock_t start = 0, end = 0;
	
	vector<vector<double>> A =
	{ {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
	{1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
	{1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
	{1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
	{1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
	{1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
	{1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
	{1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
	{1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
	{1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
	{1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
	{1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
	{1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
	{1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
	{1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
	{1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
	{1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
	{1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
	{1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
	{1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
	{1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
	{1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
	{1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
	{1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
	vector<double> b =
	{ 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
	28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
	30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
	37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 }, d = b;

	int m = b.size();
	start = clock();
	qr_decomposition(A, b, d);
	back_subs(A, b);
	end = clock();
	cout << "x0到x11的值分别为" << endl;
	for (int i = 0; i < 12; i++) {
		printf("%f\t\t", b[i]);
		if ((i + 1) % 4 == 0) {
			cout << endl;
		}
	}
	cout << endl;
	cout << "求解最小二乘运行时间为" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	for (int i = 12; i < m; i++) {
		r = r + b[i] * b[i];
	}
	r = sqrt(r);
	cout << "残向量的2范数为" << r << endl << endl;

}

void test() {
	double k = 1;
	vector<vector<double>>A = {
		{k / sqrt(2) + 1 / 2,-k / sqrt(6) + sqrt(3) / 2},
		{-k / sqrt(2),k / sqrt(6)},
		{k / sqrt(2) - 1 / 2,-k / sqrt(6) - sqrt(3) / 2} };
	vector<double> b = { 3,2,-1 }, d = b;
	qr_decomposition(A, b, d);
	back_subs(A, b);
	matrix_print(A);
	solution_print(b);
	
	

}
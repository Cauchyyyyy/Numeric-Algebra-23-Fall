#include "Exercise.h"
#include"Function.h"
#include<time.h>


void exercise_1()
{
	int N = 84; 
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	vector<int> u(N);
	vector<int> v(N);
	double max;
	clock_t start, end;      
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
	start = clock();           
	gauss_elim(A);
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();        
	solution_print(b);
	cout << "time of Gauss_elim " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;  
			

	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout <<"error "<<max<< endl;
	cout << endl;

	
	
	for (int i = 0; i < N ; i++)
	{	b[i] = 15;
		for (int j = 0; j < N ; j++) {
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
	gauss_elim_full_pivoting(A, u, v);
	vector_pb(u, b);           
	forward_subs1(A, b);
	back_subs(A, b);
	vector_qb(v, b);
	end = clock();
	solution_print(b);
	cout << "time of gauss elim full pivoting" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	
	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout << "error " << max << endl;
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
	gauss_elim_col_pivoting(A, u);
	vector_pb(u, b);	
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "time of gauss elim col pivoting" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	
	max = abs(b[0] - 1);
	for (int i = 1; i < N; i++) {
		if (max < abs(b[i] - 1))
			max = abs(1 - b[i]);
	}
	cout << "max error " << max << endl;
}

void exercise_2_1()
{
	int N = 100; 
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	clock_t start, end;

	
	for (int i = 0; i < N ; i++)
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

	srand((unsigned)time(NULL));  
	for (int i = 0; i < N; i++) {
		b[i] = rand();
	}
	cout << "original b" << endl;
	solution_print(b);

	start = clock();
	cholesky_decomp(A);
	forward_subs(A, b);
	matrix_transpose(A);   //L^T
	back_subs(A, b);
	end = clock();
	cout << "solution" << endl;
	solution_print(b);
	cout << "time of cholesky decomp" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << endl;

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

	srand((unsigned)time(NULL));  
	for (int i = 0; i < N; i++) {
		b[i] = rand();
	}
	start = clock();
	modified_cholesky_decomp(A);
	forward_subs1(A, b);
	matrix_DLT(A);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "time of modified cholesky " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
}

void exercise_2_2()
{
	int N = 40; 
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N); 
	clock_t start, end;

	

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
	cholesky_decomp(A);
	forward_subs(A, b);
	matrix_transpose(A);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "time of cholesky " << double(end - start) / CLOCKS_PER_SEC << "s" << endl<<endl;
	
	

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
	solution_print(b);
	cout << "time of modified cholesky" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
}

void exercise_3_1()
{
	int N = 100; 
	vector<vector<double>> A(N, vector<double>(N));
	vector<int> u(N);
	vector<int> v(N);
	vector<double> b(N),t(N);       
	clock_t start, end;

	
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
	
	srand((unsigned)time(NULL));  
	for (int i = 0; i < N; i++) {
		b[i] = t[i]= rand();
	}
	cout << "original b" << endl;
	solution_print(b);

	

	start = clock();
	gauss_elim(A);
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << " time of gauss elim" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << endl;

	
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
	for (int i = 0; i < N; i++) {
		b[i] = t[i];
	}

	start = clock();
	gauss_elim_full_pivoting(A, u, v);
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	vector_qb(v, b);
	end = clock();
	solution_print(b);
	cout << "time of gauss elim full pivoting" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << endl;
	
	

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
	for (int i = 0; i < N; i++) {
		b[i] = t[i];
	}

	start = clock();
	gauss_elim_col_pivoting(A, u);
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "time of gauss elim col pivoting" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
}

void exercise_3_2()
{
	int N = 40; 
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	vector<int> u(N);
	vector<int> v(N);
	clock_t start, end;

	

	for (int i = 0; i < N; i++){	
		for (int j = 0; j < N; j++){
			A[i][j] = 1 / (double(i + j + 1));
		}
	}
	for (int i = 0; i < N; i++) {
		b[i] = 0;
		for (int j = 0; j < N; j++)
		{
			b[i] = b[i] + 1 / (double(i + j + 1));
		}
	}

	start = clock();
	gauss_elim(A);
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "time of gauss elim" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << endl;

	

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			A[i][j] = 1 / (double(i + j + 1));
		}
	}
	for (int i = 0; i < N; i++) {
		b[i] = 0;
		for (int j = 0; j < N; j++)
		{
			b[i] = b[i] + 1 / (double(i + j + 1));
		}
	}

	start = clock();
	gauss_elim_full_pivoting(A,u,v);
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	vector_qb(v, b);
	end = clock();
	solution_print(b);
	cout << "time of gauss elim full pivoting" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	cout << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			A[i][j] = 1 / (double(i + j + 1));
		}
	}
	for (int i = 0; i < N; i++) {
		b[i] = 0;
		for (int j = 0; j < N; j++)
		{
			b[i] = b[i] + 1 / (double(i + j + 1));
		}
	}
	
	start = clock();
	gauss_elim_col_pivoting(A, u);		
	vector_pb(u, b);	
	forward_subs1(A, b);
	back_subs(A, b);
	end = clock();
	solution_print(b);
	cout << "time of gauss elim col pivoting" << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
}
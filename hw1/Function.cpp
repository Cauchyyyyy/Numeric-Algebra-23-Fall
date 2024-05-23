#include "Function.h"
#include <cmath>


void forward_subs(vector<vector<double>>& L, vector<double>& b)
{
	int i,j,n;
	n = b.size();
	for (i = 0; i < n - 1; i++) {
		b[i] = b[i] / L[i][i];
		for (j = i + 1; j < n; j++) {
			b[j] =b[j]- b[i] * L[j][i];
		}
	}
	b[n - 1] = b[n - 1] / L[n - 1][n - 1];
}

void forward_subs1(vector<vector<double>>& L, vector<double>& b)
{
	int i, j, n;
	n = b.size();
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			b[j]  =b[j] - b[i] * L[j][i];
		}
	}
}
	

void back_subs1(vector<vector<double>>& U, vector<double>& y)
{
	int i, j, n;
	n = y.size();
	for (i = n - 1; i > 0; i--) {
		for (j = 0; j < i; j++) {
			y[j] = y[j] - y[i] * U[j][i];
		}
	}
}

void back_subs(vector<vector<double>>& U, vector<double>& y) 
{
	int i, j, n;
	n = y.size();	
	for (i = n - 1; i > 0; i--) {
		y[i] = y[i] / U[i][i];
		for (j = 0; j < i; j++) {
			y[j] = y[j] - y[i] * U[j][i];
		}
	}
	y[0] = y[0] / U[0][0];
}

void gauss_elim(vector<vector<double>>& A)
{
	int k, i,j,n = A.size();
	for (k = 0; k < n - 1; k++) {
		for (i = k + 1; i < n; i++)
		{
			A[i][k] = A[i][k] / A[k][k];
			for (j = k + 1; j < n; j++) {
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
			}
		}
	}
}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{
	int k, i, j, m, p=0, q=0,
		n = A.size();
	double max;
	vector<double> t(n);        
	
	for (k = 0; k < n - 1; k++) {
		max = abs(A[k][k]);		
		p = q = k;         
		u[k] =v[k] = 0;           
		for (i = k; i < n; i++) {         
			for (j = k; j < n; j++) {
				if (max < abs(A[i][j])) {
					u[k] = i+1;
					p = i;
					v[k] = j+1;
					q = j;
					max = abs(A[i][j]);
				}
			}
		}
		if (max == 0) {
			printf("The matrix is strange");
			break;
		}
		else {
			if (k != p) {
				for (m = 0; m < n; m++) {           
					t[m] = A[p][m];
					A[p][m] = A[k][m];
					A[k][m] = t[m];
				}
			}
			if (k != q) {
				for (m = 0; m < n; m++)           
				{
					t[m] = A[m][q];
					A[m][q] = A[m][k];
					A[m][k] = t[m];
				}
			}
			for (i = k + 1; i < n; i++)            
			{
				A[i][k] = A[i][k] / A[k][k];
				for (j = k + 1; j < n; j++) {
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
	}
}

void gauss_elim_col_pivoting(vector<vector<double> >& A, vector<int>& u)
{
	int k, i, j, m, p = 0,s = 0, 
		n = A.size();
	double max;
	vector<double> t(n);
	
	for (k = 0; k < n - 1; k++) {
		max = abs(A[k][k]);
		p = k;
		u[k] =  0;           
		for (i = k+1; i < n; i++) {         
			if (max < abs(A[i][k])) {
					u[k] = i+1;
					p = i;
					max = abs(A[i][k]);
			}
		}		
		if (max == 0) {
			printf("The matrix is strange");
			break;
		}
		else {
			if (k != p) {
				for (m = 0; m < n; m++) {           
					t[m] = A[p][m];
					A[p][m] = A[k][m];
					A[k][m] = t[m];
				}
			}
			for (i = k + 1; i < n; i++)            
			{
				A[i][k] = A[i][k] / A[k][k];
				for (j = k + 1; j < n; j++) {
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
	}
}

void vector_pb(vector<int>& u, vector<double>& b)
{
	int i, j, n = b.size();
	double t;
	for (i = 0; i < n - 1; i++) {
		if (u[i] != 0) {
			j = u[i] - 1;
			t = b[i];
			b[i] = b[j];
			b[j] = t;
		}
	}
}

void vector_qb(vector<int>& v, vector<double>& b)
{
	int i, j, n = b.size();
	double t;
	for (i = n - 2; i >= 0; i--) {
		if (v[i] != 0) {
			j = v[i] - 1;
			t = b[i];
			b[i] = b[j];
			b[j] = t;
		}
	}
}

void cholesky_decomp(vector<vector<double>>& A)
{
	int k,i,j,n=A.size();
	for (k = 0; k < n; k++) {
		A[k][k] = sqrt(abs(A[k][k]));
		for (i = k + 1; i < n; i++) {
			A[i][k] = A[i][k] / A[k][k];
		}
		for (i = k + 1; i < n; i++) {
			for (j = i; j < n; j++)
				A[j][i] = A[j][i] - A[j][k] * A[i][k];
		}
	}
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{
	int j, i, k,
		n=A.size();
	vector<double> v(n);
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			v[i] = A[j][i] * A[i][i];
		}
		for (i = 0; i < j; i++) {
			A[j][j] = A[j][j] - A[j][i] * v[i];
		}
		for (i = j + 1; i < n; i++) {
			for (k = 0; k < j; k++) {
				A[i][j] = A[i][j] - A[i][k] * v[k];
			}
			A[i][j] = A[i][j] / A[j][j];
		}
	}
}

void matrix_DLT(vector<vector<double>>& A)
{
	matrix_transpose(A);
	int i, j, n=A.size();
	for (i = 0; i < n; i++) {
		for (j = i+1; j < n; j++) {
			A[i][j] = A[i][j] * A[i][i];
		}
	}

}

void matrix_transpose(vector<vector<double>>& A) {
	int i, j, n=A.size();
	for (i = 0; i < n; i++) {
		for (j = i+1 ; j < n; j++) {
			A[i][j] = A[j][i];
		}
	}
}

void solution_print(vector<double>& b) {
	int N = b.size();
	for (int i = 0; i < N; i++) {
		printf("%f\t\t", b[i]);
		if ((i + 1) % 4 == 0) {
			cout << endl;
		}
	}
	cout << endl;
}

void matrix_print(vector<vector<double>>& A) {
	int N = A.size();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%f\t", A[i][j]);
		}cout << endl;
	}
	cout << endl;
}
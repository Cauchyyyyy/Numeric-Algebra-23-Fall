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
	vector<double> t(n);        //交换用辅助向量
	
	for (k = 0; k < n - 1; k++) {
		max = abs(A[k][k]);		
		p = q = k;         //防止交换时因上次循环的p,q出错
		u[k] =v[k] = 0;           //对u,v清零，防止干扰
		for (i = k; i < n; i++) {         //选主元
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
				for (m = 0; m < n; m++) {           //行交换
					t[m] = A[p][m];
					A[p][m] = A[k][m];
					A[k][m] = t[m];
				}
			}
			if (k != q) {
				for (m = 0; m < n; m++)           //列交换
				{
					t[m] = A[m][q];
					A[m][q] = A[m][k];
					A[m][k] = t[m];
				}
			}
			for (i = k + 1; i < n; i++)            //Gauss变换
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
		u[k] =  0;           //对u清零，防止干扰
		for (i = k+1; i < n; i++) {         //选主元
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
				for (m = 0; m < n; m++) {           //行交换
					t[m] = A[p][m];
					A[p][m] = A[k][m];
					A[k][m] = t[m];
				}
			}
			for (i = k + 1; i < n; i++)            //Gauss变换
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
		A[k][k] = sqrt(A[k][k]);
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

void matrix_vector_multiply(vector<vector<double>>& A, vector<double>& b) {
	int n = b.size();
	vector<double>t(n);//暂时存储运算结果
	for (int i = 0; i < n; i++) {
		t[i] = 0;
		for (int j = 0 ; j < n; j++) {
			t[i] = b[j] * A[i][j] + t[i];
		}
	}
	for (int i = 0; i < n; i++) {
		b[i] = t[i];
	}
}

void matrix_full_transpose(vector<vector<double>>& A) {
	int i, j, n = A.size();
	double t;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			t = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = t;
		}
	}
}

double vectors_multiply(vector<double>& a, vector<double>& b) {
	int n=a.size();
	double t=0;
	for (int i = 0; i < n; i++) {
		t = a[i] * b[i] + t;
	}
	return t;
}

double vector_infinity_norm(vector<double>& a) {
	double max=abs(a[0]);
	int n = a.size();
	for (int i = 1; i < n; i++) {
		if (max < abs(a[i])) {
			max = abs(a[i]);
		}
	}
	return max;
}

double vector_1_norm(vector<double>& a) {
	int n = a.size();
	double t=0;
	for (int i = 0; i < n; i++) {
		t = t + abs(a[i]);
	}
	return t;
}

double matrix_inverse_infinity_norm(vector<vector<double>>& A) {
	int k = 1, n = A.size();
	double final, t1, t2;
	vector<double> x(n), w(n), v(n), z(n);
		
	for (int i = 0; i < n; i++) {
		w[i] = x[i] = static_cast<double>(1) / n;
	}
	while (k == 1) {
		matrix_full_transpose(A);
		//计算w此时的值为Bx,B=A^(-T)，解方程A^T*w=x
		forward_subs(A, w);
		back_subs1(A, w);
		for (int i = 0; i < n; i++) {
			if (w[i] == 0)
				z[i] = v[i] = 0;
			else
				z[i] = v[i] = w[i] / abs(w[i]);
		}                              //计算v=sign(w)
		matrix_full_transpose(A);
		forward_subs1(A, z);
		back_subs(A, z);//z=B^T*v
		t1 = vectors_multiply(z, x);//t1=z^T*x
		t2 = vector_infinity_norm(z);
		if (t1 >= t2) {
			final = vector_1_norm(w);
			k = 0;
		}
		else {
			for (int i = 0; i < n; i++) {
				if (t2 == abs(z[i]))
					w[i] = x[i] = 1;
				else
					w[i] = x[i] = 0;
			}
			k = 1;
		}
	}
	return final;
}

double matrix_1_norm(vector<vector<double>>& B) {
	int k = 1, n = B.size();
	double final = 0, t1 = 0, t2 = 0;
	vector<double> x(n), w(n), v(n), z(n);
	for (int i = 0; i < n; i++) {
		w[i] = x[i] = static_cast<double>(1) / n;
	}
	while (k == 1) {
		matrix_vector_multiply(B, w);//w此时的值为Bx
		
		for (int i = 0; i < n; i++) {
			if (w[i] == 0) {
				z[i] = 0;
				v[i] = 0;
			}
			else
				z[i]=v[i] = w[i] / abs(w[i]);
		}                              //计算v=sign(w)
		matrix_full_transpose(B);
		matrix_vector_multiply(B, z);//z=B^T*v
		
		t1 = vectors_multiply(z, x);//t1=z^T*x
		t2 = vector_infinity_norm(z);
		if (t1 >= t2) {
			final = vector_1_norm(w);
			k = 0;
			matrix_full_transpose(B);
		}
		else {
			for (int i = 0; i < n; i++) {
				if (t2 == abs(z[i]))
					w[i] = x[i] = 1;
				else
					w[i] = x[i] = 0;
			}
			k = 1;
			matrix_full_transpose(B);//复原矩阵B
		}
	}
	return final;
}

double matrix_1_norm_traditional(vector<vector<double>>& B) {
	double t = 0.0;
	int n = B.size();
	matrix_full_transpose(B);
	t = vector_1_norm(B[0]);
	for (int i = 1; i < n; i++) {
		if (t < vector_1_norm(B[i]))
			t = vector_1_norm(B[i]);
	}
	matrix_full_transpose(B);//复原矩阵
	return t;
}

double matrix_infinity_norm(vector<vector<double>>& A) {
	double k = 0;
	matrix_full_transpose(A);
	k = matrix_1_norm_traditional(A);
	matrix_full_transpose(A);//复原
	return k;
}
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
	n = U[0].size();
	if (n > U.size()) {
		n = U.size();
	}
	for (i = n - 1; i > 0; i--) {
		if (U[i][i] == 0) {
			cout << "the matrix is strange,U[i][i]=0,i=" <<i<< endl;
		}
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
		printf("%f\t", b[i]);
		if ((i + 1) % 6 == 0) {cout << endl;}
	}
	cout << endl<<endl;
}

void matrix_print(vector<vector<double>>& A) {
	int m = A.size(),N=A[0].size();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < N; j++) {
			printf("%f\t\t", A[i][j]);
		}cout << endl<<endl;
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
	b = t;
}

void matrix_full_transpose(vector<vector<double>>& A) {
	int i, j, m = A.size();
	double k = 0;
	
	for (i = 0; i < m; i++) {
		for (j = i + 1; j < m; j++) {
			k = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = k;
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

double householder(vector<double>&x) {
	int j = 0, n = x.size();
	double sigma = 0, beta = 0, alpha = 0,t=0;
	double k = vector_infinity_norm(x);
	
	for (j = 0; j < n; j++) {
		x[j] = x[j] / k;
	}
	vector<double> v = x;
	for (j = 1; j < n; j++) {
		sigma += x[j] * x[j];
	}
	
	if (sigma != 0) {
		alpha = sqrt(x[0] * x[0] + sigma);
		if (x[0] <= 0) {
			v[0] = x[0] - alpha;
		}
		else {
			v[0] = -sigma / (x[0] + alpha);
		}
		
		beta = 2 * v[0] * v[0] / (sigma + v[0] * v[0]);
		t = v[0];
		for (j = 0; j < n; j++) {
			x[j]=v[j] = v[j] / t;
		}
	}
	
	
	return beta;
}

void qr_decomposition(vector<vector<double>>& A,vector<double>& b, vector<double>& d) {
	int m = A.size(), n = A[0].size();
		
	for (int j = 0; j < n; j++) {
		if (j < m-1){
			vector<double> v(m-j);
			for (int i = j; i < m; i++) {
				v[i-j]=A[i][j];
			}
			
			d[j] = householder(v);
			vector<vector<double>>H(m - j, vector<double>(m - j));
			for (int i = 0; i < m - j ; i++) {
				for (int k = 0; k < m - j ; k++) {
					if (i == k) {
						H[i][k] = 1 - d[j] * v[i] * v[k];
					}
					else {
						H[i][k] =  - d[j] * v[i] * v[k];
					}
				}
			}
			vector<vector<double>>HA(m - j , vector<double>(n - j));
			vector<double> Hb(m-j);
			for (int i = 0; i < m-j; i++)
			{
				for (int k = 0; k < n - j; k++){
					HA[i][k] = A[i + j][k + j];
				}
			}
			for (int k = 0; k < m - j; k++){
				Hb[k] = b[k + j];
			}
			
			matrix_multiply(H, HA);
			matrix_vector_multiply(H, Hb);

			for (int i = 0; i < m - j; i++)
			{
				for (int k = 0; k < n - j; k++){
					 A[i + j][k + j]=HA[i][k] ;
				}
			}
			for (int k = 0; k < m - j; k++){
				 b[k + j]=Hb[k] ;
			}
			for (int i = j+1; i < m; i++){
				A[i][j] = v[i-j];
			}
			
		}
	}
}

void matrix_multiply(vector<vector<double>>& A, vector<vector<double>>& B) {
	int m = B.size(), n = B[0].size();
	vector<vector<double>>t(m, vector<double>(n));
	for (int i = 0; i < m; i++)	{
		for (int  j = 0; j < n; j++){
			t[i][j] = 0;
			for (int k = 0; k < m; k++) {
				t[i][j] = t[i][j] + A[i][k] * B[k][j];
			}
		}		
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			B[i][j] = t[i][j];
		}
	}
}

void matrix_vector_times(vector<vector<double>>& A, vector<double>& b, vector<double>& Ab) {
	int m = A.size(), n = A[0].size();
	for (int i = 0; i < m; i++) {
		Ab[i] = 0;
		for (int j = 0; j < n; j++) {
			Ab[i] = Ab[i] + A[i][j] * b[j];
		}
	}
}

void vector_vvT(vector<double>& b, vector<double>& bt, vector<vector<double>>& A) {
	int m = b.size(), n = bt.size();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = b[i] * bt[j];
		}
	}

}

void matrix_mn_transpose(vector<vector<double>>& A, vector<vector<double>>& AT) {
	int m = A.size(), n = A[0].size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			AT[i][j] = A[j][i];
		}
	}
}

void matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B) {
	int m = A.size(), n = A[0].size();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j]=A[i][j]-B[i][j];
		}
	}
}

void vector_minus(vector<double>& b, vector<double>& c) {
	int m = b.size();
	for (int i = 0; i < m; i++) {
		b[i] = b[i] - c[i];
	}
}

void qr_decomposition2(vector<vector<double>>& A, vector<double>& b, vector<double>& d) {
	int m = A.size(), n = A[0].size();

	for (int j = 0; j < n; j++) {
		if (j < m - 1) {
			vector<double> v(m - j);
			for (int i = j; i < m; i++) {
				v[i - j] = A[i][j];
			}

			d[j] = householder(v);
			
			vector<vector<double>>HA(m - j, vector<double>(n - j));
			vector<vector<double>>B(m - j, vector<double>(n - j));
			vector<vector<double>>HAT(n - j, vector<double>(m - j));
			vector<double> Hb(m - j), w(n - j);
			for (int i = 0; i < m - j; i++)
			{
				for (int k = 0; k < n - j; k++) {
					HA[i][k] = A[i + j][k + j];
				}
			}
			matrix_mn_transpose(HA, HAT);
			
			matrix_vector_times(HAT, v, w);
			for (int k = 0; k < n - j; k++) {
				w[k] = d[j] * w[k];
			}
			vector_vvT(v, w, B);
			matrix_minus(HA, B);
			for (int i = 0; i < m - j; i++)
			{
				for (int k = 0; k < n - j; k++) {
					A[i + j][k + j] = HA[i][k];
				}
			}
			for (int i = j + 1; i < m; i++) {
				A[i][j] = v[i - j];
			}
			for (int k = 0; k < m - j; k++) {
				Hb[k] = b[k + j];
			}
			double beta = vectors_multiply(v, Hb);
			for (int k = 0; k < m - j; k++) {
				v[k] = v[k] * beta * d[j];
			}
			vector_minus(Hb, v);

			
			for (int k = 0; k < m - j; k++) {
				b[k + j] = Hb[k];
			}
			

		}
	}
}

double vector_2_norm(vector<double> &a) {
	double t = 0;
	for (int i = 0; i < a.size(); i++) {
		t = t + a[i] * a[i];
	}
	t = sqrt(t);
	return t;
}

double matrix_min(vector<vector<double>>& A) {
	double t = A[0][0];
	int n = A.size();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			if (t > A[i][j]) {
				t = A[i][j];
			}
		}
	}
	return t;
}

int conjugate_gradient(vector<vector<double>>& A, vector<double>& b) {
	int n = A.size(), k = 0;
	vector<double> x(n), r(n), p(n), t(n), w(n);
	for (int i = 0; i < n; i++)
	{
		x[i] = 1;
	}
	matrix_vector_times(A, x, t);
	for (int i = 0; i < n; i++)
	{
		r[i] = b[i] - t[i];
		t[i] = 1;
	}
	double alpha=0, beta=0, tho = vectors_multiply(r, r), thoq = tho;
	while (vector_2_norm(t) >1e-7) {
		k++;
		t = x;
		if (k == 1) {
			p = r;
		}
		else {
			beta = tho / thoq;
			for (int i = 0; i < n; i++)
			{
				p[i] = r[i] + beta * p[i];
			}
		}
		matrix_vector_times(A, p, w);
		alpha = tho / vectors_multiply(p, w);
		for (int i = 0; i < n; i++)
		{
			x[i] = x[i] + alpha * p[i];
			r[i] = r[i] - alpha * w[i];
			t[i] = t[i] - x[i];
		}
		thoq = tho;
		tho = vectors_multiply(r, r);
	}
	for (int i = 0; i < n; i++)
	{
		b[i] = x[i];
	}
	return k;
}

int jacobi(vector<vector<double>>& A, vector<double>& b) {
	int k = 0, n = b.size();
	vector<vector<double>> B(n, vector<double>(n)),D(n, vector<double>(n));
	vector<double> x(n), t(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				B[i][j] = 0;
				D[i][j] = A[i][j];
			}
			else {
				B[i][j] = -A[i][j];
				D[i][j] = 0;
			}
		}
	}
	//B=(L+U)
	for (int i = 0; i < n; i++) {
		x[i] = t[i] = 1;
	}
	while (vector_2_norm(t) > 1e-7) {
		k++;
		t = x;
		matrix_vector_multiply(B, x);
		for (int i = 0; i < n; i++) {
			x[i] = (x[i] + b[i])/A[i][i];
			t[i] = t[i] - x[i];
		}
	}
	b = x;
	return k;
	
}

int g_s(vector<vector<double>>&A, vector<double>&b) {
	int k = 0, n = b.size();
	vector<vector<double>> B(n, vector<double>(n)), U(n, vector<double>(n));
	vector<double> x(n), t(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i >= j) {
				B[i][j] = A[i][j];
				U[i][j] = 0;
			}
			else {
				B[i][j] = 0;
				U[i][j] = -A[i][j];
			}
		}
	}
	//B=D-L
	forward_subs(B, b);//b=(D-L)^(-1)*b
	for (int i = 0; i < n; i++) {
		x[i] = t[i] = 1;
	}
	while (vector_2_norm(t) > 1e-7) {
		k++;
		t = x;
		matrix_vector_multiply(U, x);
		forward_subs(B, x);
		for (int i = 0; i < n; i++) {
			x[i] = x[i] + b[i];
			t[i] = t[i] - x[i];
		}
	}
	b = x;
	return k;
	
}

int sor(vector<vector<double>>& A, vector<double>& b, double omega) {
	int k = 0, n = b.size();
	vector<vector<double>> B(n, vector<double>(n)), U(n, vector<double>(n));
	vector<double> x(n), t(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i > j) {
				B[i][j] = omega * A[i][j];
				U[i][j] = 0;
			}
			else if(i==j){
				B[i][j] = A[i][j];
				U[i][j] = (1 - omega) * A[i][j];
			}
			else {
				U[i][j] = -omega * A[i][j];
				B[i][j] = 0;
			}
		}
	}
	//B=D-omega*L
	//U=(1-omega)*D+omega*U
	forward_subs(B, b);//b=(D-omega*L)^(-1)*b
	for (int i = 0; i < n; i++) {
		x[i] = t[i] = 1;
		b[i] = b[i] * omega;
	}
	while (vector_2_norm(t) > 1e-6) {
		k++;
		t = x;
		matrix_vector_multiply(U, x);
		forward_subs(B, x);
		for (int i = 0; i < n; i++) {
			x[i] = x[i] + b[i];
			t[i] = t[i] - x[i];
		}
	}
	b = x;
	return k;

}

double power_method(vector<vector<double>>& A, vector<double>& u) {
	int n = A.size();
	double t = 0, s = 0;
	int k = 0;
	
	while (1)
	{
		matrix_vector_multiply(A, u);
		t = 0;
		for (int i = 0; i < n; i++)
		{
			if (abs(u[i]) >abs(t)) {
				t = u[i];
			}
		}
		if (abs(t - s) < 1e-6) {
			break;
		}
		for (int i = 0; i < n; i++)
		{
			u[i] = u[i] / t;
		}
		s = t;
		k++;
	}
	cout << "迭代次数为" << k << endl;
	return t;
}

double max_abs_root_monic_polynomial(vector<double>& a, vector<double>& u) {
	int n = a.size();
	double k = 0;
	vector<vector<double>> A(n, vector<double>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == 0) {
				A[i][j] = - a[n - 1 - j];
			}
			else if (i == j + 1) {
				A[i][j] = 1;
			}
			else {
				A[i][j] = 0;
			}
		}
	}
	//cout << "该首一多项式的友方阵为" << endl;
	//matrix_print(A);
	return k = power_method(A,u);
}

void qr_decomposition3(vector<vector<double>>& A, vector<vector<double>>& Q) {
	int m = A.size(), n = A[0].size();
	vector<double> d(n);
	vector<vector<double>> Q1(m, vector<double>(m));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			if (i == j) {
				Q[i][j] = 1;
			}
			else {
				Q[i][j] = 0;
			}
		}
	}
	Q1 = Q;
	for (int j = 0; j < n; j++) {
		if (j < m - 1) {
			vector<double> v(m - j);
			for (int i = j; i < m; i++) {
				v[i - j] = A[i][j];
			}
			for (int i = 0; i < m; i++) {
				for (int t = 0; t < m; t++) {
					if (i == t) {
						Q1[i][t] = 1;
					}
					else {
						Q1[i][t] = 0;
					}
				}
			}
			d[j] = householder(v);
			vector<vector<double>>H(m - j, vector<double>(m - j));
			for (int i = 0; i < m - j; i++) {
				for (int k = 0; k < m - j; k++) {
					if (i == k) {
						Q1[i + j][k + j] = H[i][k] = 1 - d[j] * v[i] * v[k];
					}
					else {
						Q1[i + j][k + j] = H[i][k] = -d[j] * v[i] * v[k];
					}
				}
			}
			vector<vector<double>>HA(m - j, vector<double>(n - j));
			
			for (int i = 0; i < m - j; i++)
			{
				for (int k = 0; k < n - j; k++) {
					HA[i][k] = A[i + j][k + j];
				}
			}
			matrix_multiply(H, HA);
			matrix_multiply(Q1, Q);
			for (int i = 0; i < m - j; i++)
			{
				for (int k = 0; k < n - j; k++) {
					A[i + j][k + j] = HA[i][k];
				}
			}
			
		}
	}
	matrix_full_transpose(Q);
}

void direct_qr_algorithm(vector<vector<double>>& A,vector<vector<double>> &Q,int k) {
	int n = A.size();
	for (int i = 0; i < k; i++) {
		qr_decomposition3(A, Q);
		matrix_multiply2(A, Q);
		
	}
	
}

void matrix_multiply2(vector<vector<double>>& A, vector<vector<double>>& B) {
	int m = A.size(), n = A[0].size();
	vector<vector<double>>t(m, vector<double>(n));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			t[i][j] = 0;
			for (int k = 0; k < n; k++) {
				t[i][j] = t[i][j] + A[i][k] * B[k][j];
			}
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = t[i][j];
		}
	}
}

void up_hessenberg_decomposition(vector<vector<double>>& A, vector<vector<double>>& Q) {
	int n = A.size();
	vector<double> d(n);
	vector<vector<double>> Q1(n, vector<double>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				Q[i][j] = 1;
			}
			else {
				Q[i][j] = 0;
			}
		}
	}
	Q1 = Q;
	for (int k = 0; k < n-2; k++) {
		
		vector<double> v(n - k - 1);
		for (int i = k+1; i < n; i++) {
				v[i - k - 1] = A[i][k];
		}
		for (int i = 0; i < n; i++) {
			for (int t = 0; t < n; t++) {
				if (i == t) {
						Q1[i][t] = 1;
					}
				else {
						Q1[i][t] = 0;
					}
				}
			}
		d[k] = householder(v);
		vector<vector<double>>H(n - k - 1, vector<double>(n - k - 1));
		for (int i = 0; i < n-k-1; i++) {
			for (int j = 0; j < n-k-1; j++) {
				if (i == j) {
						Q1[i + k+1][k+1 + j] = H[i][j] = 1 - d[k] * v[i] * v[j];
					}
				else {
						Q1[i + k+1][k+1 + j] = H[i][j] = -d[k] * v[i] * v[j];
					}
				}
			}
		vector<vector<double>>HA(n - k - 1, vector<double>(n - k));

		for (int i = 0; i < n-k-1; i++)
			{
			for (int j = 0; j < n-k; j++) {
					HA[i][j] = A[i + k + 1][j + k];
				}
			}
		matrix_multiply(H, HA);
		matrix_multiply(Q1, Q);
		for (int i = 0; i < n - k - 1; i++){
			for (int j = 0; j < n - k; j++) {
					A[i + k + 1][k + j] = HA[i][j];
			}
		}
		vector<vector<double>>AH(n, vector<double>(n - k-1));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - k - 1; j++) {
					AH[i][j] = A[i][j + k + 1];
			}
		}
		matrix_multiply2(AH, H);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - k - 1; j++) {
				A[i][j + k + 1]=AH[i][j];
			}
		}
	}
}

void francis_displacement(vector<vector<double>>& A, vector<vector<double>>& Q) {
	int a = A.size(), n = a - 1, m = n - 1;
	vector<vector<double>> Q1(a, vector<double>(a));
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < a; j++) {
			Q[i][j] = 0;
		}Q[i][i] = 1;
	}
	Q1 = Q;
	if (a <= 2) {
		direct_qr_algorithm(A, Q,1);
	}
	else {
		double s = A[m][m] + A[n][n],
			t = A[m][m] * A[n][n] - A[m][n] * A[n][m],
			x = A[0][0] * A[0][0] + A[0][1] * A[1][0] - s * A[0][0] + t,
			y = A[1][0] * (A[0][0] + A[1][1] - s),
			z = A[1][0] * A[2][1];
		
		for (int k = -1; k <=n - 3; k++) {
			for (int i = 0; i < a; i++) {
				for (int j = 0; j < a; j++) {
					Q1[i][j] = 0;
				}Q1[i][i] = 1;
			}
			vector<double> v(3);
			v = { x,y,z };
			double beta = householder(v);
			int q = max(0, k);
			vector<vector<double>> H(3, vector<double>(3));
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					Q1[i + k + 1][j + k + 1] = H[i][j] = -beta * v[i] * v[j];
				}
				H[i][i] += 1;
				Q1[i + k + 1][i + k + 1] = H[i][i];
			}
			matrix_multiply(Q1, Q);
			vector<vector<double>> HA(3, vector<double>(n-q+1));
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < a - q; j++) {
					HA[i][j] = A[k + 1 + i][j + q];
				}
			}
			matrix_multiply(H, HA);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < a - q; j++) {
					 A[k + 1 + i][j + q]=HA[i][j];
				}
			}
			int r = min(k + 4, n);
			vector<vector<double>> AH(r + 1, vector<double>(3));
			for (int i = 0; i < r+1; i++) {
				for (int j = 0; j < 3; j++) {
					AH[i][j] = A[i][k + 1 + j];
				}
			}
			matrix_multiply2(AH,H);
			for (int i = 0; i < r + 1; i++) {
				for (int j = 0; j < 3; j++) {
					 A[i][k + 1 + j]=AH[i][j];
				}
			}
			x = A[k + 2][k + 1];
			y = A[k + 3][k + 1];
			if (k < n - 3) {
				z = A[k + 4][k + 1];
			}
		}
		for (int i = 0; i < a; i++) {
			for (int j = 0; j < a; j++) {
				Q1[i][j] = 0;
			}
			Q1[i][i] = 1;
		}
		vector<double> v(2);
		v = { x,y };
		double beta = householder(v);
		vector<vector<double>> H = 
		{ {1 - beta * v[0] * v[0],-beta * v[0] * v[1]},
			{-beta * v[1] * v[0],1 - beta * v[1] * v[1]} };

		Q1[m][m] = H[0][0]; Q1[m][n] = H[0][1]; Q1[n][m] = H[1][0]; Q1[n][n] = H[1][1];
		matrix_multiply(Q1, Q);
		vector<vector<double>> HA(2, vector<double>(3));
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 3; j++) {
				HA[i][j] = A[a - 2 + i][a - 3 + j];
			}
		}
		matrix_multiply(H, HA);
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 3; j++) {
				A[a - 2 + i][a - 3 + j]=HA[i][j];
			}
		}
		vector<vector<double>> AH(a, vector<double>(2));
		for (int i = 0; i < a; i++) {
			for (int j = 0; j < 2; j++) {
				AH[i][j] = A[i][a - 2 + j];
			}
		}
		matrix_multiply2(AH, H);
		for (int i = 0; i < a; i++) {
			for (int j = 0; j < 2; j++) {
				 A[i][a - 2 + j]=AH[i][j];
			}
		}
	}
}

int implicit_qr_algorithm(vector<vector<double>>& A, vector<vector<double>>& Q) {
	int n = A.size(), m = 0;
	int k = 0;
	vector<vector<double>>Q1(n, vector<double>(n));
	up_hessenberg_decomposition(A, Q);
	while (m < n) {
		k++;//if (k > n*n) { break; }
		int l = 0, t = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				Q1[i][j] = 0;
			}
			Q1[i][i] = 1;
		}
		for (int i = 1; i < n; i++) {
			if (abs(A[i][i - 1]) < 1e-10) {
				A[i][i - 1] = 0;
			}
		}
		//matrix_print(A);
	label:
		for (int i = n - m - 1; i > 0; i--) {
			if (A[1][0] == 0) {//下 次对角线全为0
				m = n;
			}
			if (A[i][i - 1] != 0) {
				m = n - i - 1;
				t = i;
				break;
			}
		}
		if (m == n||m==n-1) {
			break;
		}
		if (t == 1) {
			if (pow(A[0][0] + A[1][1], 2)
				- 4 * (A[0][0] * A[1][1] - A[0][1] * A[1][0]) < 0) {
				break;
			}
		}
		else {
			if (A[t - 1][t - 2] == 0) {//处理2阶阵
				if (pow(A[t][t] + A[t - 1][t - 1], 2)
					- 4 * (A[t][t] * A[t - 1][t - 1] - A[t][t - 1] * A[t - 1][t]) < 0) {
					m += 2;
					goto label;//2阶方阵，特征值是共轭复根
				}
				else {
					l = t - 1;
				}

			}
			else {//找首个为0的次对角线元
				for (int i = t - 2; i > 0; i--) {
					if (A[i][i - 1] == 0) {
						l = i; break;
					}
				}
			}
		}
		//处理H22
		int a = n - l - m;
		vector<vector<double>>H22(n - l - m, vector<double>(n - l - m));
		
		vector<vector<double>>H23(n - l - m, vector<double>(m));
		vector<vector<double>>P(n - l - m, vector<double>(n - l - m));
		for (int i = 0; i < a; i++) {
			for (int j = 0; j < a; j++) {
				H22[i][j] = A[l + i][l + j];
			}
		}
		
		for (int i = 0; i < a; i++) {
			for (int j = 0; j < m; j++) {
				H23[i][j] = A[l + i][n - m + j];
			}
		}
		francis_displacement(H22, P);

		for (int i = 0; i < a; i++) {
			for (int j = 0; j < a; j++) {
				A[l + i][l + j] = H22[i][j];
			}
		}
		for (int i = 0; i < a; i++) {
			for (int j = 0; j < a; j++) {
				Q1[l + i][l + j] = P[i][j];
			}
		}
		matrix_multiply(P, H23);  
		for (int i = 0; i < a; i++) {
			for (int j = 0; j < m; j++) {
				A[l + i][n - m + j] = H23[i][j];
			}
		}
		if (l != 0) { 
			vector<vector<double>>H12(l, vector<double>(n - l - m));
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < a; j++) {
					H12[i][j] = A[i][l + j];
				}
			}
			matrix_full_transpose(P);
			matrix_multiply2(H12, P);
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < a; j++) {
					A[i][l + j] = H12[i][j];
				}
			}
		}
		
		matrix_multiply(Q1, Q);
		//matrix_multiply(Q1, A); matrix_full_transpose(Q1); matrix_multiply2(A, Q1);
	}return k;
}

void find_eigenvalue_from_schur(vector<vector<double>>& A) {
	int n = A.size();
	for (int i = n - 1; i > 0; i--) {
		if (i != 0) {
			if (abs(A[i][i - 1]) <1e-10) {
				cout << A[i][i]<<endl;
			}
			else {
				double b = -A[i][i] - A[i - 1][i - 1],
					c = A[i][i] * A[i - 1][i - 1] - A[i - 1][i] * A[i][i - 1];
				double im = sqrt(4 * c - b * b) / 2, re = -b / 2;
				cout << re << "+" << im << "i"<<endl;
				cout << re << "-" << im << "i"<<endl;
				i--;
			}
		}
	}if (abs(A[1][0]) < 1e-10) {
		cout << A[0][0] << endl;
	}
}

double bisection(vector<vector<double>>& A, int k) {
	int n = A.size();
	vector<double>x(n), y(n);
	for (int i = 0; i < n; i++) {
		x[i] = A[i][i];
		if (i > 0) {
			y[i] = A[i][i - 1];
		}
	}
	y[0] = 0;
	double u = matrix_infinity_norm(A) + 1, l = -u, miu = 0;
	int t = 0;
	
	while (abs(l - u) > 1e-10) {
		t++;
		int s = 0; 
		double q = x[0] - miu;
		for (int i = 0; i < n; i++) {
			if (q < 0) {
				s++;
			}
			if (i < n - 1) {
				if (q == 0) {
					q = 1e-10;
				}
				q = x[i + 1] - miu - y[i + 1] * y[i + 1] / q;
			}
		}
		if (s >= k) {
			u = miu;
		}
		else {
			l = miu;
		}
		miu = (l + u) / 2;
	}
	cout << "二分法迭代次数为" << t << endl;
	return miu;
}

void inverse_power_method(vector<vector<double>>& A, double t, vector<double>& u,int l) {
	int n = A.size();
	double t1=0,s = 0;
	int k = 0;
	vector<vector<double>> A1(n, vector<double>(n));
	vector<int> v(n);
	vector<double>x(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A1[i][j] = A[i][j];
		}A1[i][i] = A[i][i] - t;
	}
	gauss_elim_col_pivoting(A1, v);
	while (k<l)
	{
		k++;
		vector_pb(v, u);
		forward_subs1(A1, u);
		back_subs(A1, u);

		t1 = vector_2_norm(u);
		
		
		for (int i = 0; i < n; i++)
		{
			u[i] = u[i] / t1;
			x[i] = x[i] - u[i];
		}
		if (abs(t1 - s) < 1e-6|| vector_2_norm(x)<1e-6) {
			break;
		}
		s = t1;
		
		
	}
	cout << "迭代次数为" << k << endl;
	
}

void inverse_power_method2(vector<vector<double>>& A, double t, vector<double>& u) {
	int n = A.size();
	double t1 = 0, s = 0;
	int k = 0;
	vector<vector<double>> A1(n, vector<double>(n));
	vector<int> v(n);
	vector<double>x(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A1[i][j] = A[i][j];
		}A1[i][i] = A[i][i] - t;
	}
	gauss_elim_col_pivoting(A1, v);
	while (1)
	{
		k++;
		vector_pb(v, u);
		forward_subs1(A1, u);
		back_subs(A1, u);

		t1 = vector_2_norm(u);

		
		for (int i = 0; i < n; i++)
		{
			u[i] = u[i] / t1;
			x[i] = abs(u[i])-abs(x[i]);
		}
		if (abs(t1 - s) < 1e-6|| vector_2_norm(x) < 1e-6) {
			break;
		}
		s = t1;


	}
	cout << "迭代次数为" << k << endl;

}

void jacobi_method(vector<vector<double>>& A) {
	int n = A.size();
	double delta = 1;
	
	int k = 0, a = 0;
	vector<vector<double>>B = A;
	while (delta > 1e-10) {
		a++;
		k = 0;
		for (int p = 0; p < n - 1; p++) {
			for (int q = p + 1; q < n; q++) {
				if (abs(A[p][q]) > delta) {
					k++;
					double tao = 0, t = 0, s = 0, c = 0, sgn = 0;
					tao = (A[q][q] - A[p][p]) / (2 * A[p][q]);
					if (tao > 0) {
						sgn = 1;
					}
					if (tao < 0) {
						sgn = -1;
					}
					t = sgn / (abs(tao) + sqrt(1 + tao * tao));
					if (tao == 0) {
						t = 1;
					}
					c = 1 / sqrt(1 + t * t);
					s = t * c;
					for (int i = 0; i < n; i++) {
						if (i != p && i != q) {
							A[i][p] = A[p][i] = c * B[i][p] - s * B[i][q];
							A[i][q] = A[q][i] = s * B[i][p] + c * B[i][q];
						}
					}
					A[p][p] = c * c * B[p][p] - 2 * s * c * B[p][q] + s * s * B[q][q];
					A[q][q] = c * c * B[q][q] + 2 * s * c * B[p][q] + s * s * B[p][p];
					A[p][q] = A[q][p] = (c * c - s * s) * B[p][q] + s * c * (B[p][p] - B[q][q]);
					B = A;
				}
			}
		}
		if (k == 0) {
			delta = delta / n;
		}
		//cout << delta << endl;matrix_print(A);
	}
	cout << "扫描次数为" << a << endl;
}

void sort_s_to_l(vector<double>& a) {
	int n = a.size();
	double t = 0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n - 1-j; i++) {
			if (a[i] > a[i + 1]) {
				t = a[i];
				a[i] = a[i + 1];
				a[i + 1] = t;
			}
		}
	}
}

void jacobi_method2(vector<vector<double>>& A, vector<vector<double>>& Q) {
	int n = A.size();
	double delta = 1;

	int k = 0, a = 0;
	vector<vector<double>>B = A, Q1 = Q;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Q[i][j] = 0;
		}Q[i][i] = 1;
	}Q1 = Q;
	while (delta > 1e-10) {
		a++;
		k = 0;
		for (int p = 0; p < n - 1; p++) {
			for (int q = p + 1; q < n; q++) {
				if (abs(A[p][q]) > delta) {
					
					k++;
					double tao = 0, t = 0, s = 0, c = 0, sgn = 1;
					tao = (A[q][q] - A[p][p]) / (2 * A[p][q]);
					
					if (tao < 0) {
						sgn = -1;
					}
					t = sgn / (abs(tao) + sqrt(1 + tao * tao));
					
					c = 1 / sqrt(1 + t * t);
					s = t * c;
					
					for (int i = 0; i < n; i++) {
						if (i != p && i != q) {
							A[i][p] = A[p][i] = c * B[i][p] - s * B[i][q];
							A[i][q] = A[q][i] = s * B[i][p] + c * B[i][q];
						}
						Q[i][p] = Q1[i][p] * c - Q1[i][q] * s;
						Q[i][q] = Q1[i][q] * c + Q1[i][p] * s;
					}
					A[p][p] = c * c * B[p][p] - 2 * s * c * B[p][q] + s * s * B[q][q];
					A[q][q] = c * c * B[q][q] + 2 * s * c * B[p][q] + s * s * B[p][p];
					A[p][q] = A[q][p] = (c * c - s * s) * B[p][q] + s * c * (B[p][p] - B[q][q]);
					
					B = A;
					Q1 = Q;
				}
			}
		}
		if (k == 0) {
			delta = delta / n;
		}
		//cout << delta << endl;matrix_print(A);
	}
	//cout << "Jacobi法扫描次数为" << a << endl;
}

void bidiagonalisation(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V) {
	int m = A.size(), n = A[0].size();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			U[i][j] = 0;
		}U[i][i] = 1;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			V[i][j] = 0;
		}V[i][i] = 1;
	}
	vector<vector<double>> U1(m, vector<double>(m));
	vector<vector<double>> V1(n, vector<double>(n));
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				U1[i][j] = 0;
			}U1[i][i] = 1;
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				V1[i][j] = 0;
			}V1[i][i] = 1;
		}
		vector<double> v(m - k);
		for (int i = 0; i < m - k; i++) {
			v[i] = A[i + k][k];
		}
		double beta = householder(v);
		for (int i = 0; i < m-k; i++) {
			for (int j = 0; j < m-k; j++) {
				U1[i + k][j + k] = -beta * v[i] * v[j];
			}U1[i + k][i + k] += 1;
		}
		matrix_multiply(U1, A);
		matrix_multiply(U1, U);
		if (k < n - 2) {
			vector<double> u(n - k - 1);
			for (int i = 0; i < n - k - 1; i++) {
				u[i] = A[k][i + 1 + k];
			}
			beta = householder(u);
			for (int i = 0; i < n-k-1; i++) {
				for (int j = 0; j < n-k-1; j++) {
					V1[i + k + 1][j + k + 1] = -beta * u[i] * u[j];
				}V1[i + k + 1][i + k + 1] += 1;
			}
			matrix_multiply2(A,V1);
			matrix_multiply2(V,V1);
		}
	}
}

void wilkinson_displacement(vector<vector<double>>&A, vector<vector<double>>& U, vector<vector<double>>& V) {
	int n = A.size();
	if (n <= 2) {
		//由第七章习题13
		vector<vector<double>> C(2, vector<double>(2));
		double t = ( A[1][0] - A[0][1])/(A[0][0] + A[1][1]);
		C[0][0] = C[1][1] = 1 / sqrt(t * t + 1);
		C[0][1] = t / sqrt(1 + t * t);
		C[1][0] = -C[0][1];
		matrix_multiply(C, A);
		jacobi_method2(A, V);
		//matrix_full_transpose(V);
		C[1][0] = -C[1][0];	C[0][1] = -C[0][1];
		U = V;
		matrix_multiply(C, U);
		
		
	}
	else {
		vector<double> d(n), g(n - 1);
		
		for (int i = 0; i < n - 1; i++) {
			d[i] = A[i][i];
			g[i] = A[i][i + 1];
		}d[n - 1] = A[n - 1][n - 1];
		double alpha = d[n - 1] * d[n - 1] + g[n - 2] * g[n - 2],
			delta = (d[n - 2] * d[n - 2] + g[n - 3] * g[n - 3] - alpha) / 2,
			beta = d[n - 2] * g[n - 2],
			sign = 0;
		if (delta > 0) { sign = 1; }
		else if (delta < 0) { sign = -1; }
		double miu = alpha - beta * beta / (delta + sign * sqrt(beta * beta + delta * delta)),
			y = d[0] * d[0] - miu, z = d[0] * g[0];
		double c = 1, s = 0, t = 0, temp = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				U[i][j] = 0;
				V[i][j] = 0;
			}U[i][i] = 1;
			V[i][i] = 1;
		}
		
		for (int k = 0; k < n - 1; k++) {
			
			if (y != 0) {
				t = -z / y;
				s = t / sqrt(t * t + 1), c = 1 / sqrt(t * t + 1);
			}
			else {
				s = 1, c = 0;
			}
			
			for (int i = 0; i < n; i++) {
				temp = V[i][k];
				V[i][k] = c * V[i][k] - s * V[i][k + 1];
				V[i][k + 1] = s * temp + c * V[i][k + 1];
			}
			
			if (k > 0) {
				g[k - 1] = c * y - s * z;
			}
			y = c * d[k] - s * g[k];
			z = -s * d[k + 1];
			g[k] = s * d[k] + c * g[k];
			d[k+1] = c * d[k + 1];
			if (y != 0) {
				t = -z / y;
				s = t / sqrt(t * t + 1), c = 1 / sqrt(t * t + 1);
			}
			else {
				s = 1, c = 0;
			}
			for (int i = 0; i < n; i++) {
				temp = U[i][k];
				U[i][k] = c * U[i][k] - s * U[i][k + 1];
				U[i][k + 1] = s * temp + c * U[i][k + 1];
			}
			
			
			d[k] = c * y - s * z;
			if (k < n - 2) {
				y = c * g[k] - s * d[k + 1];
				z = -s * g[k + 1];
				d[k + 1] = s * g[k] + c * d[k + 1];
				g[k + 1] = c * g[k + 1];
			}
			else {
				temp = g[k];
				g[k] = c * g[k] - s * d[k + 1];
				d[k + 1] = s * temp + c * d[k + 1];
			}
			
			

		}
		for (int i = 0; i < n - 1; i++) {
			 A[i][i]=d[i];
			 A[i][i + 1]=g[i];
		} A[n - 1][n - 1]=d[n-1];
	}
}

void svd_decomposition(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V) {
	int m = A.size(), n = A[0].size();
	int pp = 0;
	bidiagonalisation(A, U, V);
	double epsilon = 1e-10;
	matrix_full_transpose(U);
	while (pp<300) {
		pp++;
		vector<vector<double>> U1(m, vector<double>(m));
		vector<vector<double>> V1(n, vector<double>(n));
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				U1[i][j] = 0;
			}U1[i][i] = 1;
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				V1[i][j] = 0;
			} V1[i][i] = 1;
		}
		label:
		for (int i = 0; i < m; i++) {
			for (int j = i; j < i+2&&j<n; j++) {
				if (abs(A[i][j]) < epsilon) {
					A[i][j] = 0;
				}
			}
		}
		int p = 0, q = n, t = 0;
		for (int i = n - 2; i >= 0; i--) {
			if (A[i][i + 1] != 0) {
				t = i;
				q = n - i - 2;
				break;
			}
		}
		if (q == n || q == n - 1) {
			break;
		}
		for (int i = t - 1; i >= 0; i--) {
			if (A[i][i + 1] == 0) {
				p = i + 1;
				break;
			}
		}
		
		if (n - p - q <=1) {
			cout << "something wrong" << endl;
			break;
		}
		t = -1;
		vector<vector<double>> B(n - p - q, vector<double>(n - p - q));
		for (int i = 0; i < n - p - q; i++) {
			for (int j = 0; j < n - p - q; j++) {
				B[i][j] = A[i + p][j + p];
			}
			if (B[i][i] == 0 && i != n - p - q - 1) {
				t = i;
			}
		}
		if (t < 0 ) {
			vector<vector<double>> P(n - p - q, vector<double>(n - p - q));
			vector<vector<double>> Q(n - p - q, vector<double>(n - p - q));
			
			//matrix_print(B);
			wilkinson_displacement(B, P, Q);
			for (int i = 0; i < n - p - q; i++) {
				for (int j = 0; j < n - p - q; j++) {
					A[i + p][j + p] = B[i][j];
					U1[i + p][j + p] = P[i][j];
					V1[i + p][j + p] = Q[i][j];
				}
			}
			matrix_multiply2(U, U1);
			matrix_multiply2(V, V1);
		}
		else {
			cout << "still need coding" << endl;
			break;
			goto label;
		}
	}
	for (int i = 0; i < n; i++) {
		if (A[i][i] < 0) {
			for (int j = 0; j < n; j++) {
				V[j][i] = -V[j][i];
				A[j][i] = -A[j][i];
			}
		}
	}
	cout << "SVD迭代次数为" << pp << endl;
}
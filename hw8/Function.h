#pragma once
#include<iostream>
#include<vector>
using namespace std;

void matrix_transpose(vector<vector<double>>& A);//transpose of triangular matrix

void forward_subs(vector<vector<double>>& L, vector<double>& b);//forward substitution

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//forward substitution for diagonal elements are 1

void back_subs(vector<vector<double>>& U, vector<double>& b);//back substitution

void back_subs1(vector<vector<double>>& U, vector<double>& y);//back substitution for diagonal elements are 1

void gauss_elim(vector<vector<double>>& A);//Gauss elimination

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//full pivoting Gauss elimination

void gauss_elim_col_pivoting( vector<vector<double>>& A, vector<int>& u);//column pivoting Gauss elimination

void vector_pb(vector<int>&u,vector<double>&b);//compute vector P*b

void vector_qb(vector<int>& v, vector<double>& b);//compute vector Q*b

void cholesky_decomp(vector<vector<double>>& A);//Cholesky decomposition

void modified_cholesky_decomp(vector<vector<double>>& A);//modified Cholesky decomposition

void matrix_DLT(vector<vector<double>>& A);//compute D*L^T

void solution_print(vector<double>& b);//print the solution

void matrix_print(vector<vector<double>>& A);//print the matrix

void matrix_vector_multiply(vector<vector<double>>& A, vector<double>& b);
//compute the product of matrix A and vector b
//store the result in b

void matrix_full_transpose(vector<vector<double>>& A);//transpose

double vectors_multiply(vector<double>& a, vector<double>& b);
//compute the product of two vectors
//return a*.b (inner product a.k.a dot product)

double vector_infinity_norm(vector<double>& a);//infinity norm of vector a

double vector_1_norm(vector<double>& a);//1 norm of vector a

double vector_2_norm(vector<double>& a);//2 norm of vector a

double matrix_1_norm(vector<vector<double>>& B);//�Ż������ƾ����1����

double matrix_1_norm_traditional(vector<vector<double>>& B);//һ�㷽����������1����

double matrix_infinity_norm(vector<vector<double>>& A);//infinity norm of matrix A

double matrix_inverse_infinity_norm(vector<vector<double>>& A);//���ƾ���A^(-1)�������

double householder(vector<double>&x);//��xͨ��Householder�任Ϊe_1����ʽ,���ئ�=2/(v^T*v)

void qr_decomposition(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//����A��QR�ֽ�,������ΪA��������ΪR����������ÿ��householder�任v��ֵ��d�洢v[0]

void matrix_multiply(vector<vector<double>>& A, vector<vector<double>>& B);
//A \in R^{m*m}, B \in R^{m*n}, B = A*B

void matrix_vector_times(vector<vector<double>>& A, vector<double>& b, vector<double>& Ab);
//A \in R^{m*n}, b \in R^{n*1}, Ab \in R^{m*1}

void matrix_mn_transpose(vector<vector<double>>& A, vector<vector<double>>& AT);
//AT=A^T

void qr_decomposition2(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//����A��QR�ֽ�,������ΪA��������ΪR����������ÿ��householder�任v��ֵ��d�洢v[0]
//�����þ���˻�

void matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B);
//A=A-B

void vector_vvT(vector<double>& b, vector<double>& bt, vector<vector<double>>& A);
//A=b*b^T

void vector_minus(vector<double>& b, vector<double>& c);
//b=b-c

double matrix_min(vector<vector<double>>& A);
//minimum element of matrix A

int conjugate_gradient(vector<vector<double>>& A, vector<double>& b);
//ʵ�ù����ݶȷ����Ax=b,��洢��b,���ص�������

int jacobi(vector<vector<double>>& A, vector<double>& b);
//Jacobi iterative method, return the number of iteration times
//solution stored in b

int g_s(vector<vector<double>>& A, vector<double>& b);
//Gauss-Seidel iterative method, return the number of iteration times

int sor(vector<vector<double>>& A, vector<double>& b, double omega);
//Successive Over Relaxation iterative method, return the number of iteration times

double power_method(vector<vector<double>>& A, vector<double>& u);
//�ݷ�������ֵ���ģ,u�ǵ����ĳ�ʼֵ

double max_abs_root_monic_polynomial(vector<double>& a, vector<double>& u);
//�����ϵ��Ϊa����������һ����ʽ��ģ����

void qr_decomposition3(vector<vector<double>>& A, vector<vector<double>>& Q);
//QR�ֽ⣬A=QR��R�洢��A��

void direct_qr_algorithm(vector<vector<double>>& A, vector<vector<double>>& Q, int k);
//ֱ�ӵ�QR�㷨

void up_hessenberg_decomposition(vector<vector<double>>& A, vector<vector<double>>& Q);
//�������Hessenberg�ֽ⣬Q�Ǳ��εı任����,�õ����¾���QAQ^T�洢��A

int implicit_qr_algorithm(vector<vector<double>>& A, vector<vector<double>>& Q);
//��ʽQR�㷨

void matrix_multiply2(vector<vector<double>>& A, vector<vector<double>>& B);
//�������˻�A*B������A��m*n�ģ�B��n*n�ģ�����洢��A

void francis_displacement(vector<vector<double>>& A, vector<vector<double>>& Q);
//�Ծ���A��˫�ز�λ�Ƶ�QR�����㷨������洢��A��Q�Ǳ������ڱ任��������

void find_eigenvalue_from_schur(vector<vector<double>>& A);
//��ʵSchur��׼��������������ֵ

double bisection(vector<vector<double>>& A, int k);
//���ַ���Գ����ԽǾ���Ĵ�С�����k������ֵ

void inverse_power_method(vector<vector<double>>& A, double t, vector<double>& u,int l);
//���ݷ���A������ֵt��Ӧ����������

void inverse_power_method2(vector<vector<double>>& A, double t, vector<double>& u);

void jacobi_method(vector<vector<double>>& A);
//����Jacobi����

void jacobi_method2(vector<vector<double>>& A, vector<vector<double>>& Q);
//����Jacobi�����������Ӧ����������

void sort_s_to_l(vector<double>& a);
//������a��С��������

void bidiagonalisation(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V);
//����Ķ��Խǻ�

void wilkinson_displacement(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V);
//Wilkinsonλ�Ƶ�SVD����

void svd_decomposition(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V);
//��ʽQR����ԭ����е�����ֵ�ֽ�

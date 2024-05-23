#pragma once
#include<iostream>
#include<vector>
using namespace std;

void matrix_transpose(vector<vector<double>>& A);//���Ǿ���ת��

void forward_subs(vector<vector<double>>& L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>>& U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>>& U, vector<double>& y);//�Խ�ԪΪ1�Ļش���

void gauss_elim(vector<vector<double>>& A);//Gauss��ȥ��

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//ȫ��ԪGauss��ȥ��

void gauss_elim_col_pivoting( vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��

void vector_pb(vector<int>&u,vector<double>&b);//��������P*b����ѡ��

void vector_qb(vector<int>& v, vector<double>& b);//��������Q*b����ѡ��

void cholesky_decomp(vector<vector<double>>& A);//�Գ��������׼Cholesky�ֽ�

void modified_cholesky_decomp(vector<vector<double>>& A);//�Ľ���ƽ������

void matrix_DLT(vector<vector<double>>& A);//�������D*L^T����ѡ��

void solution_print(vector<double>& b);//�����

void matrix_print(vector<vector<double>>& A);//�������

void matrix_vector_multiply(vector<vector<double>>& A, vector<double>& b);
//��������A*b,�洢��b

void matrix_full_transpose(vector<vector<double>>& A);//һ�㷽��ת��

double vectors_multiply(vector<double>& a, vector<double>& b);//����a,b�ڻ�

double vector_infinity_norm(vector<double>& a);//����a�������

double vector_1_norm(vector<double>& a);//����a��1����

double vector_2_norm(vector<double>& a);//����a��2����

double matrix_1_norm(vector<vector<double>>& B);//�Ż������ƾ����1����

double matrix_1_norm_traditional(vector<vector<double>>& B);//һ�㷽����������1����

double matrix_infinity_norm(vector<vector<double>>& A);//�������������

double matrix_inverse_infinity_norm(vector<vector<double>>& A);//���ƾ���A^(-1)�������

double householder(vector<double>&x);//��xͨ��Householder�任Ϊe_1����ʽ,���ئ�=2/(v^T*v)

void qr_decomposition(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//����A��QR�ֽ�,������ΪA��������ΪR����������ÿ��householder�任v��ֵ��d�洢v[0]

void matrix_multiply(vector<vector<double>>& A, vector<vector<double>>& B);
//�������˻�A*B������A��m*m�ģ�B��m*n�ģ�����洢��B

void matrix_vector_times(vector<vector<double>>& A, vector<double>& b, vector<double>& Ab);
//����A*b��A��m*n�ģ�b��n*1�ģ��洢��Ab��

void matrix_mn_transpose(vector<vector<double>>& A, vector<vector<double>>& AT);
//һ������ת�ã��洢��AT��

void qr_decomposition2(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//����A��QR�ֽ�,������ΪA��������ΪR����������ÿ��householder�任v��ֵ��d�洢v[0]
//�����þ���˻�

void matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B);
//A-B�洢��A

void vector_vvT(vector<double>& b, vector<double>& bt, vector<vector<double>>& A);
//v*v^T�洢��A

void vector_minus(vector<double>& b, vector<double>& c);
//b=b-c

double matrix_min(vector<vector<double>>& A);
//����Ԫ�ص���С����

int conjugate_gradient(vector<vector<double>>& A, vector<double>& b);
//ʵ�ù����ݶȷ����Ax=b,��洢��b,���ص�������

int jacobi(vector<vector<double>>& A, vector<double>& b);
//Jacobi����������洢��b�����ص�������

int g_s(vector<vector<double>>& A, vector<double>& b);
//G-S������

int sor(vector<vector<double>>& A, vector<double>& b, double omega);
//SOR������

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
//�Ծ���A��˫�ز�λ���㷨������洢��A��Q�Ǳ������ڱ任��������

void find_eigenvalue_from_schur(vector<vector<double>>& A);
//��ʵSchur��׼��������������ֵ
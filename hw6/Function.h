#pragma once
#include<iostream>
#include<vector>
using namespace std;

void matrix_transpose(vector<vector<double>>& A);//三角矩阵转置

void forward_subs(vector<vector<double>>& L, vector<double>& b);//前代法

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//对角元为1的前代法

void back_subs(vector<vector<double>>& U, vector<double>& b);//回代法

void back_subs1(vector<vector<double>>& U, vector<double>& y);//对角元为1的回代法

void gauss_elim(vector<vector<double>>& A);//Gauss消去法

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//全主元Gauss消去法

void gauss_elim_col_pivoting( vector<vector<double>>& A, vector<int>& u);//列主元Gauss消去法

void vector_pb(vector<int>&u,vector<double>&b);//计算向量P*b【可选】

void vector_qb(vector<int>& v, vector<double>& b);//计算向量Q*b【可选】

void cholesky_decomp(vector<vector<double>>& A);//对称正定阵标准Cholesky分解

void modified_cholesky_decomp(vector<vector<double>>& A);//改进的平方根法

void matrix_DLT(vector<vector<double>>& A);//计算矩阵D*L^T【可选】

void solution_print(vector<double>& b);//输出解

void matrix_print(vector<vector<double>>& A);//输出矩阵

void matrix_vector_multiply(vector<vector<double>>& A, vector<double>& b);
//计算向量A*b,存储于b

void matrix_full_transpose(vector<vector<double>>& A);//一般方阵转置

double vectors_multiply(vector<double>& a, vector<double>& b);//向量a,b内积

double vector_infinity_norm(vector<double>& a);//向量a的无穷范数

double vector_1_norm(vector<double>& a);//向量a的1范数

double vector_2_norm(vector<double>& a);//向量a的2范数

double matrix_1_norm(vector<vector<double>>& B);//优化法估计矩阵的1范数

double matrix_1_norm_traditional(vector<vector<double>>& B);//一般方法计算矩阵的1范数

double matrix_infinity_norm(vector<vector<double>>& A);//计算矩阵的无穷范数

double matrix_inverse_infinity_norm(vector<vector<double>>& A);//估计矩阵A^(-1)的无穷范数

double householder(vector<double>&x);//将x通过Householder变换为e_1的形式,返回β=2/(v^T*v)

void qr_decomposition(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//矩阵A的QR分解,计算结果为A的上三角为R，下三角是每次householder变换v的值，d存储v[0]

void matrix_multiply(vector<vector<double>>& A, vector<vector<double>>& B);
//计算矩阵乘积A*B，其中A是m*m的，B是m*n的，结果存储于B

void matrix_vector_times(vector<vector<double>>& A, vector<double>& b, vector<double>& Ab);
//计算A*b，A是m*n的，b是n*1的，存储于Ab中

void matrix_mn_transpose(vector<vector<double>>& A, vector<vector<double>>& AT);
//一般矩阵的转置，存储于AT中

void qr_decomposition2(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//矩阵A的QR分解,计算结果为A的上三角为R，下三角是每次householder变换v的值，d存储v[0]
//不采用矩阵乘积

void matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B);
//A-B存储于A

void vector_vvT(vector<double>& b, vector<double>& bt, vector<vector<double>>& A);
//v*v^T存储于A

void vector_minus(vector<double>& b, vector<double>& c);
//b=b-c

double matrix_min(vector<vector<double>>& A);
//矩阵元素的最小分量

int conjugate_gradient(vector<vector<double>>& A, vector<double>& b);
//实用共轭梯度法求解Ax=b,解存储于b,返回迭代次数

int jacobi(vector<vector<double>>& A, vector<double>& b);
//Jacobi迭代法，解存储于b，返回迭代次数

int g_s(vector<vector<double>>& A, vector<double>& b);
//G-S迭代法

int sor(vector<vector<double>>& A, vector<double>& b, double omega);
//SOR迭代法

double power_method(vector<vector<double>>& A, vector<double>& u);
//幂法求特征值最大模,u是迭代的初始值

double max_abs_root_monic_polynomial(vector<double>& a, vector<double>& u);
//求各项系数为a各分量的首一多项式的模最大根

void qr_decomposition3(vector<vector<double>>& A, vector<vector<double>>& Q);
//QR分解，A=QR，R存储于A中

void direct_qr_algorithm(vector<vector<double>>& A, vector<vector<double>>& Q, int k);
//直接的QR算法

void up_hessenberg_decomposition(vector<vector<double>>& A, vector<vector<double>>& Q);
//方阵的上Hessenberg分解，Q是本次的变换矩阵,得到的新矩阵QAQ^T存储于A

int implicit_qr_algorithm(vector<vector<double>>& A, vector<vector<double>>& Q);
//隐式QR算法

void matrix_multiply2(vector<vector<double>>& A, vector<vector<double>>& B);
//计算矩阵乘积A*B，其中A是m*n的，B是n*n的，结果存储于A

void francis_displacement(vector<vector<double>>& A, vector<vector<double>>& Q);
//对矩阵A做双重步位移算法，结果存储在A，Q是本次用于变换的正交阵

void find_eigenvalue_from_schur(vector<vector<double>>& A);
//从实Schur标准型输出矩阵的特征值
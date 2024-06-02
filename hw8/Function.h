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

double matrix_1_norm(vector<vector<double>>& B);//优化法估计矩阵的1范数

double matrix_1_norm_traditional(vector<vector<double>>& B);//一般方法计算矩阵的1范数

double matrix_infinity_norm(vector<vector<double>>& A);//infinity norm of matrix A

double matrix_inverse_infinity_norm(vector<vector<double>>& A);//估计矩阵A^(-1)的无穷范数

double householder(vector<double>&x);//将x通过Householder变换为e_1的形式,返回β=2/(v^T*v)

void qr_decomposition(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//矩阵A的QR分解,计算结果为A的上三角为R，下三角是每次householder变换v的值，d存储v[0]

void matrix_multiply(vector<vector<double>>& A, vector<vector<double>>& B);
//A \in R^{m*m}, B \in R^{m*n}, B = A*B

void matrix_vector_times(vector<vector<double>>& A, vector<double>& b, vector<double>& Ab);
//A \in R^{m*n}, b \in R^{n*1}, Ab \in R^{m*1}

void matrix_mn_transpose(vector<vector<double>>& A, vector<vector<double>>& AT);
//AT=A^T

void qr_decomposition2(vector<vector<double>>& A, vector<double>& b, vector<double>& d);
//矩阵A的QR分解,计算结果为A的上三角为R，下三角是每次householder变换v的值，d存储v[0]
//不采用矩阵乘积

void matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B);
//A=A-B

void vector_vvT(vector<double>& b, vector<double>& bt, vector<vector<double>>& A);
//A=b*b^T

void vector_minus(vector<double>& b, vector<double>& c);
//b=b-c

double matrix_min(vector<vector<double>>& A);
//minimum element of matrix A

int conjugate_gradient(vector<vector<double>>& A, vector<double>& b);
//实用共轭梯度法求解Ax=b,解存储于b,返回迭代次数

int jacobi(vector<vector<double>>& A, vector<double>& b);
//Jacobi iterative method, return the number of iteration times
//solution stored in b

int g_s(vector<vector<double>>& A, vector<double>& b);
//Gauss-Seidel iterative method, return the number of iteration times

int sor(vector<vector<double>>& A, vector<double>& b, double omega);
//Successive Over Relaxation iterative method, return the number of iteration times

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
//对矩阵A做双重步位移的QR迭代算法，结果存储在A，Q是本次用于变换的正交阵

void find_eigenvalue_from_schur(vector<vector<double>>& A);
//从实Schur标准型输出矩阵的特征值

double bisection(vector<vector<double>>& A, int k);
//二分法求对称三对角矩阵的从小到大第k个特征值

void inverse_power_method(vector<vector<double>>& A, double t, vector<double>& u,int l);
//反幂法求A的特征值t对应的特征向量

void inverse_power_method2(vector<vector<double>>& A, double t, vector<double>& u);

void jacobi_method(vector<vector<double>>& A);
//过关Jacobi方法

void jacobi_method2(vector<vector<double>>& A, vector<vector<double>>& Q);
//过关Jacobi方法，保存对应的特征向量

void sort_s_to_l(vector<double>& a);
//对序列a从小到大排序

void bidiagonalisation(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V);
//矩阵的二对角化

void wilkinson_displacement(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V);
//Wilkinson位移的SVD迭代

void svd_decomposition(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V);
//隐式QR迭代原理进行的奇异值分解

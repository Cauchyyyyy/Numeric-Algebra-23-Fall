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
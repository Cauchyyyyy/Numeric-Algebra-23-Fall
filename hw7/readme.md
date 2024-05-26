# Homework 7

## Description
This is the solutions for Homework 7 of the Numeric Algebra course.

## Installation
1. Visual Studio 2022

## Usage
open the homework7.sln and click to run

## Introduction
### 第1题
求实对称三对角阵的全部特征值和特征向量。

(1) 用 C++ 编制利用过关 Jacobi 方法求实对称三对角阵全部特征值和特征向量的通用子程序。

(2) 利用你所编制的子程序求50, 60, 70, 80, 90, 100 阶矩阵
$$
A=\begin{bmatrix}
\;4 & 1 & 0 & 0 & \cdots & 0 \;\\
\;1 & 4 & 1 & 0 & \cdots & 0 \;\\
\;0 & 1 & 4 & 1 & \cdots & 0 \;\\
\;\vdots & & \ddots & \ddots & \ddots \;\\
\;0 & \cdots & 0 & 1 & 4 & 1\;\\
\;0 & 0 & \cdots & 0 & 1 & 4\; 
\end{bmatrix}
$$
的全部特征值和特征向量。

源文件中的输出格式：

n=xx, 

Ak=

[矩阵]

Qk=

[矩阵]

迭代次数：x, 用时 xxx s.
### 第2题
求实对称三对角阵的指定特征值及对应的特征向量.

(1) 用 C++ 编制先利用二分法求实对称三对角阵指定特征值，再利用反幂法求对应特征向量的通用子程序。

(2) 利用你所编制的子程序求100阶矩阵

$$
A=\begin{bmatrix}
\;2 & -1 & 0 & 0 & \cdots & 0 \;\\
\;-1 & 2 & -1 & 0 & \cdots & 0 \;\\
\;0 & -1 & 2 & -1 & \cdots & 0 \;\\
\;\vdots & & \ddots & \ddots & \ddots \;\\
\;0 & \cdots & 0 & -1 & 2 & -1\;\\
\;0 & 0 & \cdots & 0 & -1 & 2\; 
\end{bmatrix}
$$
的最大和最小特征值及对应的特征向量
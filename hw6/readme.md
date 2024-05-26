# Homework 6

## Description
This is the solutions for Homework 6 of the Numeric Algebra course.

## Installation
1. Visual Studio 2022

## Usage
open the homework6.sln and click to run

## Introduction
### 第1题
求多项式方程的模最大根。

(1) 编制利用幂法求多项式方程 
$$f(x) = x^n + \alpha_{n-1}x^{n-1} +\cdots+ \alpha_1 x + \alpha_0 = 0 $$
的模最大根的通用子程序。

(2) 利用编制的子程序求下列各高次方程的模最大根。

(i) $x^3 + x^2 - 5x + 3 = 0;$

(ii) $x^3 - 3x - 1 = 0;$

(iii) $x^8 + 101x^7 + 208.01x^6 + 10891.01x^5 + 9802.08x^4 + 79108.9x^3-99902x^2 + 790x-1000 = 0.$


### 第2题
(1) 编制利用隐式 QR 算法 (课本算法 6.4.3) 求一个实矩阵的全部特征值的通用子程序。

(2) 利用编制的子程序计算方程 
$$x^{41} + x^3 + 1 = 0 $$

的全部根。

(3)设
$$A=\begin{bmatrix}
    \;9.1 & 3.0 & 2.6 & 4.0 \;\\
    \;4.2 & 5.3 & 4.7 & 1.6 \;\\
    \;3.2 & 1.7 & 9.4 & x \;\\
    \;6.1 & 4.9 & 3.5 & 6.2\; 
    \end{bmatrix}
$$

求当 x = 0.9, 1.0, 1.1 时 A 的全部特征值.
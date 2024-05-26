# Homework 5

## Description
This is the solutions for Homework 5 of the Numeric Algebra course.

## Installation
1. Visual Studio 2022

## Usage
open the homework5.sln and click to run

## Introduction
### 第1题
考虑 Dirichlet 问题
$$ \left\{
\begin{aligned}
& -\Delta u+u=f(x,y),(x,y)\in [0,1]×[0,1] \\
& u|_{\Gamma}=\varphi
\end{aligned}
\right.
$$

其中 $\Gamma$ 为正方形区域的边界。类似于模型问题，我们得到差分方程
$$ \left\{
\begin{aligned}
& (1 + \frac{h^2}{4})u_{i,j}-\frac{1}{4}(u_{i-1,j}+u_{i,j-1}+u_{i+1,j}+u_{i,j+1}) = \frac{h^2}{4}f_{i,j},i,j=1,\cdots ,n-1\\
& u_{i,0}=\phi_{i,0},u_{i,n}=\phi_{i,n},i=0,1,\cdots ,n \\
& u_{0,j}=\phi_{0,j},u_{n,j}=\phi_{n,j},j=0,1,\cdots ,n
\end{aligned}
\right.
$$

按照自然顺序排列得到系数矩阵为
$$A=\begin{bmatrix}
    \;S' & B &  &  &  &  \;\\
    \;B & S' & B &  &  &  \;\\
    \; & B & S' & B &  &  \;\\
    \;  & & \ddots & \ddots & \ddots \;\\
    \; &  &  & B & S' & B\;\\
    \; &  &  &  & B & S'\; 
    \end{bmatrix}
    $$



其中 $B = -I/4, I$ 为 n-1 阶单位矩阵，$S'$是对角元均为 $1+h^2/4$, 
次对角元均为-1/4 的n-1 阶对称三对角阵。对 $f(x, y) = sin(xy), \varphi(x, y) = x^2 + y^2, n = 20$。

用共轭梯度法求解差分方程，要求输出迭代次数、求解所用时间和解向量，迭代终止条件为  
$||x_{k+1}-x_k||_{\infty} < 10^{-7}$。
注意边界条件与线性方程组的关系!!

### 第2题
用 Hilbert 矩阵测试你所编写的共轭梯度法程序：
$$a_{ij}=\frac{1}{i+j-1}, b_i =\frac{1}{3}\sum_{j=1}^{n}a_{ij} , 0 \leq i, j \leq n-1$$
对 n = 20, 40, 60, 80 分别求解，观察解是否准确，迭代停止条件自定，输出迭代次数、求解所用时间和解向量。

### 第3题
分别用 Jacobi 迭代法，G-S 迭代法和共轭梯度法求解下述方程，
输出迭代次数、求解所用时间和解向量，并对结果给出解释。
$$\begin{pmatrix}
    \;10 & 1 & 2 & 3 & 4 \;\\
    \;1 & 9 & -1 & 2 & -3 \;\\
    \;2 & -1 & 7 & 3 & -5 \;\\
    \;3 & 2 & 3 & 12 & -1 \;\\
    \;4 & -3 & -5 & -1 & 15\; 
    \end{pmatrix}
    \begin{pmatrix}
        \;x_1\;\\
        \;x_2\;\\
        \;x_3\;\\
        \;x_4\;\\
        \;x_5\; 
        \end{pmatrix}
    =\begin{pmatrix}
        \;12\;\\
        \;-27\;\\
        \;14\;\\
        \;-17\;\\
        \;12\; 
        \end{pmatrix}
$$

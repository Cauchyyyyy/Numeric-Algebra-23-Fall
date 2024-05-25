# Homework 2

## Description
This is the solutions for Homework 2 of the Numeric Algebra course.

## Installation
1. Visual Studio 2022

## Usage
open the homework2.sln and click to run

## Instruction
将算法 2.5.1 编写成通用的子程序, 然后用你编写的程序完成以下计算任务:

（1）估计并输出 5 到 20 阶 Hilbert 矩阵的 ∞ 范数条件数


（2）设
$$A_n=
    \left(
    \begin{array}{cccccc}
     1 & 0 & 0 & 0 & 0 & 1 \\
     -1 & 1 & 0 & 0 & 0 & 1 \\
     -1 & -1 & 1 & 0 & 0 & 1 \\
     \vdots  & \vdots  & \ddots & \ddots & \ddots & \vdots  \\
     -1 & \ddots & \ddots & -1 & 1 & 1  \\
     -1 & -1 & \cdots  & \cdots  & -1 & 1 \\
    \end{array}
    \right)
$$
    
随机选取  $ x \in \mathbb{R}$, 并计算出$ b = A_n x$, 然后用列主元高斯消去法求解该方程组。
    估计 n 从 5 到 30 时计算解 x̂ 的精度，并与真实相对误差作比较,输出真实相对误差
    和估计的相对误差上界。

## Contributing
Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

#   Numeric Algebra Homework1

##  problem

1. 用不选主元的 guass 消去法、全主元 guass 消去法和列主元 guass 消去法求解 84 阶方程组:
$$ \left(
\begin{array}{cccccc}
 6 & 1 &   &   &   &   \\
 8 & 6 & 1 &  &   &   \\
  & 8 & 6 & 1 &  &  \\
  &  & 8 & 6 & \ddots &  \\
  &  &  & \ddots & \ddots & 1 \\
  &  &  &  & 8 & 6 \\
\end{array}
\right)\left(
\begin{array}{c}
 x_1 \\
 x_2 \\
 x_3 \\
 \vdots  \\
 \vdots  \\
 x_{84} \\
\end{array}
\right)=\left(
\begin{array}{c}
 7 \\
 15 \\
 15 \\
 \vdots  \\
 \vdots  \\
 14 \\
\end{array}
\right)
$$
输出计算结果，计算结果和准确解的误差以及运行时间。


2. 用平方根法和改进平方根法求解对称正定方程组 Ax=b:

(1)b 随机选取, A 为 100 阶矩阵
$$
\left(
\begin{array}{cccccc}
 10 & 1 &   &   &   &   \\
 1 & 10 & 1 &   &   &   \\
   & 1 & 10 & 1 &   &   \\
   &   & 1 & 10 & \ddots &   \\
   &   &   & \ddots & \ddots & 1 \\
   &   &   &   & 1 & 10 \\
\end{array}
\right)$$

(2) 矩阵A为 40 阶 Hilbert 矩阵, 其中$a_{\text{ij}}=\frac{1}{i+j-1}$，
向量 b 的第 i 个分量$b_i=\sum _{j=1}^n \frac{1}{i+j-1}$


3. 用不选主元的 guass 消去法、全主元 guass 消去法和列主元 guass 消去法求解2
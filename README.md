# HermiteInterpolationCppCode


&emsp;
## 1 Hermite插值方法简介
简单来说，给定 $M$ 个互异的点，以及这些点上的函数值以及一阶导数值，共计 $N$ 个插值条件（ $N$ 可以不等于 $2M$ ），则可以确定唯一的一个 $N-1$ 个插值多项式，称为Hermite插值多项式。
### 1.1 插值方法
首先构造一个Hermite插值多项式：
$$H\left( x \right) = {a_0} + {a_1}x + {a_{\rm{2}}}{x^{\rm{2}}}{\rm{ + }} \cdots  + {a_{n - 1}}{x^{n - 1}} = \sum\limits_{i = 0}^{n - 1} {{a_i}{x^i}}  
\tag{1}$$
$$\frac{{\partial H\left( x \right)}}{{\partial x}} = {a_1} + 2{a_{\rm{2}}}x{\rm{ + }} \cdots  + \left( {n - 1} \right){a_{n - 1}}{x^{n - 2}} = \sum\limits_{i = 1}^{n - 1} {i{a_i}{x^{i - 1}}} \tag{2} $$

假设有 $N_{1}$ 个函数值条件以及 $N_{2}$ 个一阶导数条件($N = {N_1} + {N_2}$)，将这 $N$ 个条件分别带入到式(1)以及式(2),即可求解系数 ${a_0},{a_1}, \cdots ,{a_{n - 1}}$。

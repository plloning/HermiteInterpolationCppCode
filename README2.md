&emsp;
# Hermite插值
姓名：叶科奇，学号：23111605，日期：23.10.7

&emsp;
## 1 Hermite插值方法简介
简单来说，给定$M$个互异的点，以及这些点上的函数值以及一阶导数值，共计$N$个插值条件（$N$可以不等于$2M$），则可以确定唯一的一个$N-1$个插值多项式，称为Hermite插值多项式。
### 1.1 插值方法
首先构造一个Hermite插值多项式：
$$H\left( x \right) = {a_0} + {a_1}x + {a_{\rm{2}}}{x^{\rm{2}}}{\rm{ + }} \cdots  + {a_{n - 1}}{x^{n - 1}} = \sum\limits_{i = 0}^{n - 1} {{a_i}{x^i}}  
\tag{1}$$
$$\frac{{\partial H\left( x \right)}}{{\partial x}} = {a_1} + 2{a_{\rm{2}}}x{\rm{ + }} \cdots  + \left( {n - 1} \right){a_{n - 1}}{x^{n - 2}} = \sum\limits_{i = 1}^{n - 1} {i{a_i}{x^{i - 1}}} \tag{2} $$

假设有$N_{1}$个函数值条件以及$N_{2}$个一阶导数条件（$N = {N_1} + {N_2}$），将这$N$个条件分别带入到式(1)以及式(2),即可求解系数${a_0},{a_1}, \cdots ,{a_{n - 1}}$。
&emsp;
## 2 算例
<div align=center>
<image src ='123.jpeg' >
</div>
题目中12个插值条件可以唯一确定一个11阶Hermite插值多项式(全局插值)，也可以分段，每段得到一个三阶多项式（题设要求）。分段插值的计算结果如下表

| $x$ | 0.5 | 1.5 |2.5 |3.5 |4.8 |
| :----:| :----: | :----: |:----: |:----: |:----: |
| $H(x)$ | 0.8125 | 0.3075 | 0.1375 |0.075372 |0.0415869 |
| $f(x)$ | 0.8 | 0.307692 | 0.137931 |0.0754717 |0.0415973 |
| $error$ | 0.0156 | -0.000625 | -0.003215 |-0.00132137 |-0.000251802 |

注意上表中的误差取相对误差，其定义为
$$error=\frac{H(x)-f(x)}{f(x)} $$
## 3 改进
### 3.1 方法1
在$[0,1]$之间再取一个点 $x=0.6$，构造一个5阶多项式。新的计算结果如下

| $x$ | 0.5 | 
| :----:| :----: | 
| $H(x)$ | 0.799795 | 
| $f(x)$ | 0.8 | 
| $error$ | -0.000256812 | 

如果额外点的位置为 $x=0.1$，那么

| $x$ | 0.5 | 
| :----:| :----: | 
| $H(x)$ | 0.797275 | 
| $f(x)$ | 0.8 | 
| $error$ | -0.00340653 | 

可以看到，插值点距离测量点距离越近，效果越好。
### 3.2 方法2
取切比雪夫零点作为节点，可以降低余项，避免龙格现象，这在插值点个数较多时尤为重要。
## 4 代码
代码使用了Eigen库的矩阵计算模块。
### 4.1 4阶插值函数
```cpp
void HermiteInterMetod(Eigen::MatrixXd &a, Eigen::MatrixXd &ax, Eigen::MatrixXd &x, Eigen::MatrixXd &N)
{

    Eigen::MatrixXd Q(4, 4);
    Eigen::MatrixXd invQ = Q;
    Eigen::MatrixXd F(4, 1);
    F.block(0, 0, 2, 1) = a;
    F.block(2, 0, 2, 1) = ax;

    double temp = 0.0;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            temp = pow(x(i, 0), j);
            Q(i, j) = temp;
        }
    }
    for (int i = 0; i < 2; i++)
    {
        for (int j = 1; j < 4; j++)
        {
            temp = j * pow(x(i, 0), j - 1);
            Q(i + 2, j) = temp;
        }
    }
    invQ = Q.inverse();
    N = invQ * F;
};
```
### 4.1 6阶插值函数
```cpp
void HermiteInterMetodRefine(Eigen::MatrixXd &a, Eigen::MatrixXd &ax, Eigen::MatrixXd &x, Eigen::MatrixXd &N)
{

    Eigen::MatrixXd Q(6, 6);
    Eigen::MatrixXd invQ = Q;

    // Eigen::MatrixXd f(N1, 1);
    // Eigen::MatrixXd fx(N2, 1);
    Eigen::MatrixXd F(6, 1);
    F.block(0, 0, 3, 1) = a;
    F.block(3, 0, 3, 1) = ax;

    double temp = 0.0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            temp = pow(x(i, 0), j);
            Q(i, j) = temp;
        }
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 1; j < 6; j++)
        {
            temp = j * pow(x(i, 0), j - 1);
            Q(i + 3, j) = temp;
        }
    }
    invQ = Q.inverse();
    N = invQ * F;
    cout << Q << endl;
    cout << "-----------" << endl;
    cout << invQ << endl;
    cout << "-----------" << endl;
    cout << F << endl;
    cout << "-----------" << endl;
};
```
### 4.1 主函数
```cpp
#include <iostream>
#include <Eigen/Dense>
using namespace std;

double Calf(double x)
{
    return 1 / (1 + pow(x, 2));
}
double Calfx(double x)
{
    return -2 * x / pow((1 + pow(x, 2)), 2);
}
int main()
{
    int Num = 3;
    Eigen::MatrixXd a(Num, 1);
    Eigen::MatrixXd ax(Num, 1);
    Eigen::MatrixXd x(Num, 1);
    Eigen::MatrixXd N(4, 1);
    double tempa = 0;
    double tempc = 0.1;
    double tempb = 1;
    double t = 0.5;
    a(0, 0) = Calf(tempa);
    ax(0, 0) = Calfx(tempa);
    x(0, 0) = tempa;
    a(1, 0) = Calf(tempc);
    ax(1, 0) = Calfx(tempc);
    x(1, 0) = tempc;
    a(2, 0) = Calf(tempb);
    ax(2, 0) = Calfx(tempb);
    x(2, 0) = tempb;

    Eigen::MatrixXd tp(6, 1);

    for (int i = 0; i < 6; i++)
    {
        tp(i, 0) = pow(t, i);
    }
    HermiteInterMetodRefine(a, ax, x, N);
    Eigen::MatrixXd h = tp.transpose() * N;
    cout << h(0, 0) << endl;
    cout << Calf(t) << endl;
    cout << (h(0, 0) - Calf(t)) / Calf(t) << endl;

}
```

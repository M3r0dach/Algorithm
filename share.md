<script type="text/javascript"
   src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
# 数学方面的有关算法及应用

标签（空格分隔）： 数论 离散数学

---

## 前言
&emsp;&emsp;本次分享的算法和数学知识可以这么理解，都是可以用整数类型的变量可以实现的

## 一些预备知识和符号说明
1. 取模 mod
2. 对实数x向下取整  `[x]` `c语言里正整数除法自动向下取整`
3. a在群中的逆元 $a^{-1}$
4. 位运算 与：&， 右移位:>>， 左移位:<<
5. 复杂度符号 O(f(n))
6. 大整数类型  LL  `typedef long long LL;` 或者 `typedef __int128_t LL;`  
7. 求和符号 $\Sigma$

## 数论方面

### 一个基础的算法:(快速幂/模重复平方算法)  

* 递归版思想  
:    在对p取模的情况下，计算 $a^x$的结果时不需要通过 $a^{x-1}$，而是可以直接通过 $a^{[x/2]}$ 得出  
    代码 O(log n)
```c++
LL fast_pow(LL a, LL x, LL mod) {
    if(x==1) return a;
    LL ret=fast_pow(a,x/2);
    ret = ret*ret%mod;
    if(x&1) ret=ret*a%mod;
    return ret;
}
```
* 非递归版思想
:    把x用二进制表示，拆解成若干个2的k次方向加  
则 $x= \sum\limits_{k=0}^{64} (x>>k\&1)*(1<<k)$  
则 $a^x=\prod\limits_{k=0}^{64}(x>>k\&1)*a^{1<<k}$  
    代码 复杂度O(log n)  

```c++
LL fast_pow(LL a, LL x, LL mod) {
    LL ret=1;
    while(x) {
        if(x&1) ret=ret*a%mod;
        a=a*a%mod;
        x>>=1;
    }
    return ret;
}
```
* 应用
:    求一个有理数p/q的第k位小数  
$$
\begin{eqnarray}ans
		&=&[p/q*10^k] \mod 10\\
		&=&[p*10^k/q] \mod 10
	\end{eqnarray}\\
对于[a/b] \mod{c} 有 \\
[a/b]\mod{c} = [(a \mod{bc}) /b] \mod {c} \\
所以 ans =[(p*10^k \mod{10q}) /q] \mod 10
$$
* 题目
:    已知  $\pi=\sum\limits_{k=0}^\infty{\frac{1}{{{{16}^k}}}\left({\frac{4}{{8k+1}}-\frac{2}{{8k+4}}-\frac{1}{{8k+5}}-\frac{1}{{8k+6}}}\right)}$
求 $\pi$十六进制下的小数第k位  
    + 讲讲做法,做为作业

---

### 求解素数模下的逆元
  * 在线算法
:    + 费马小定理
```c++
\\复杂度O(log p)
LL inv(LL a, LL p) {
    return fast_pow(a, p-2, p);        
}
```
    + 扩展欧几里得
  * 离线算法 求1-n关于p的逆元
:    借用扩展欧几里得算法的思想  
$$
若 xa+yp=1 有解\\
则 pm+(p\mod a)n=1有解\\
即pm+(p-[p/a]*a)n=1\\
即-[p/a]na+(n+m)p=1\\
所以可以取x=-[p/a]*inv(p\mod a)$$

```c++
LL inv[maxn];
const LL p = 1e9+7;
void init(int n) {
    inv[1]=1;
    for(int i=2; i<=n; ++i)
        inv[i] = -(p/i)*inv[p%i]%p;
}
```
  * 应用  
:    已知$a,b,p$ 求解$x,使得a^x = b \mod p$  
 两种方法  
    + 普通解法,复杂度**O(p)**
    + 第二种,**BayStep_GiantStep**  
$$
任选一个正整数m,若有解必有x=ms+r,0<=r<m\\
则a^{ms+r}=b,即a^r=b*a^{-ms}\\
生成字典[a^0:0,a^1:1,...a^{m-1}:m-1]\\
枚举s,需枚举p/m次\\
查找b*a^{-ms}是否在字典中,若在则得出答案\\
复杂度O(p/m*log(m)+m*log(m))\\
取m=\sqrt{p}有最优值O(\sqrt(p)*log(p))
$$

---

### 平方剩余
* 问题描述
:    $给定a,和奇素数p,求解方程x^2=n \pmod p$

* 判定是否有解,勒让德符号  
$$
\left(\frac{n}{p}\right)=a^{(p-1)/2}=
\begin{cases}
1,&n\text{在模$p$意义下是二次剩余}\\
-1,&n\text{在模$p$意义下是非二次剩余}\\
0,&n\equiv0\pmod p
\end{cases}
$$
* 普通解法 **O(p/2)**
* Cipolla O(log(p))
    1. 不断地随机出一个数a，使得 $\left(\frac{a^2-n}{p}\right)=-1$
    2. 给($a^2-n$)一个 $域\mathbb {F}_{p^2}$ 让它可以开根号（类比$\sqrt{−1}$所在的复数域）。我们将 $\sqrt{a^2-n}$ 定义为这个域的“虚数单位元”，设它为ω
    3. $x\equiv (a+\omega)^\frac{p+1}{2}\pmod p$, 注意是 $域\mathbb {F}_{p^2}$ 下的快速幂
    4. 证明 <http://blog.csdn.net/a_crazy_czy/article/details/51959546>

---

## 极客大挑战两道题正解讲解
### 出题的缘由和想法
### misc 找规律
* 题意
:    小C是个科研狗，一天他获得了一组数列，是
0, 1, 1, 2, 8, 18, 59, 155, 460, 1276, 3672, 10357, 29533
他觉得他能根据已知信息确定规律，比如这组数列的第30项。
那聪明的你也知道答案吗？注意项数是从第0项开始
* 思考
:    数列的类型有很多种,比如多项式函数,或者递推数列,他们的数字特征都不一样,如果推断是k次多项式最好的方法是开个k次方根,看得到的数据是否线性,反正这个数列呈指数分布,是个递推数列
* 线性递推数列的一般解法
:    设$\{a_n\}$ 是个m阶多项式
$$
a_n = a_{n-1}x_1+a_{n-2}x_2+a_{n-3}x_3....a_{n-m}x_m\\
则\\
a_{n+1} = a_{n}x_1+a_{n-1}x_2+a_{n-2}x_3....a_{n-m+1}x_m\\
a_{n+2} = a_{n+1}x_1+a_{n}x_2+a_{n-1}x_3....a_{n-m+2}x_m\\
a_{n+3} = a_{n+2}x_1+a_{n+1}x_2+a_{n}x_3....a_{n-m+3}x_m\\
...\\
a_{n+m} = a_{n+m-1}x_1+a_{n+m-2}x_2+a_{n+m-3}x_3....a_{n}x_m\\
求解方程组得出答案\\
一些小优化,由于低阶递推式在高阶有无数解\\
所以方程可能有多解,可以先取个很大的m得到一个矩阵\\
它的行列式值应该是0,在求这个矩阵的秩就OK了\\
$$
matlab代码
```matlab
a=[0, 1, 1, 2, 8, 18, 59, 155, 460, 1276, 3672, 10357, 29533]
y=a(9:13)'
A=[a(8:-1:4); a(9:-1:5); a(10:-1:6); a(11:-1:7);a(12:-1:8)]
x=rank(A)
A=A(1:x,1:x)
y=y(1:x)
inv(A)*y
```
* 应用
    :    在很多计数问题或者寻找满足某个规律的第n个数字时都非常有用
    1. 例题1
        :    海伦数a, 三角形满足三边长为 a-1,a,a+1,且三角形面积为整数,求第1000个海伦数
    2. 作业
        :    在$n*4$的矩阵里，填满$1*2$的骨牌，有多少种可能性.
        <http://acm.hdu.edu.cn/showproblem.php?pid=1992>

### code 排列
* 题意
:    求一个排列的第前k排列是什么?
```python
p=[14,1,13,2,12,15,9,8,11,6,5,10,3,7,4]`
k=307674368000;
```
* 解法康托展开
:    对于一个集合有n!个排列,按照字典序可以对这n!个排列进行排序,则每个排列和它的排名一一对应,即$f(P)=rank, P = f^{-1}(rank)$,$f 和 f^{-1}$可以通过康托展开和康托逆展开实现
* ***应用*** 八数码<http://acm.hdu.edu.cn/showproblem.php?pid=1043>

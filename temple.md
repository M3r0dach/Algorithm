<h1 align = "center">目录</h1>

1. DP
    - 斜率优化的DP
    - 插头DP
    - 数位DP
    - 大背包
1. 计算几何  
    - 凸包
    - 动态凸包
    - 旋转卡壳法
    - 矩形面积求并
    - Geo模板
1. 数论
    - 离散对数
    - 欧拉函数
    - 素数测试
    - 素数个数
    - Lucas定理
    - 平方剩余
1. 区间问题
    - 莫队算法
    - fzu2224
    - Treap
    - 树上莫队cot2
    - Splay
    - HJTree
    - huafen
    - TreeinTee
1. 离散数学
	- FFT
    - NTT
    - FWT
    - 拆系数FFT
    - 线性递推
	- 高精度
1. 图论
    - LCA tarjian
1. 字符串
    - 后缀数组
    - 最长回文子串
    - ShiftOr
    - LCP
    - 后缀自动机
    - 最小表示法


* * *

# DP
## 斜率优化的DP

```c++
#include<cstdio>
typedef long long ll;
const int MAXN = 5e5+1e2;
struct Pos{
    ll x, y;
    Pos(ll xx=0, ll yy=0)
    :x(xx), y(yy){}
    Pos operator - (const Pos& b) {
        return Pos(x-b.x, y-b.y);
    }
}stk[MAXN];
int n, m, sz;
ll sum[MAXN], dp[MAXN];
ll cross(const Pos& a, const Pos& b) {
    return a.x*b.y-a.y*b.x;
}
int main() {
    while(~scanf("%d%d", &n, &m)) {
        for(int i=1; i<=n; ++i) {
            scanf("%I64d", sum+i);
            sum[i] += sum[i-1];
        }
        sz = 0;
        stk[sz++] = Pos(0,0);
        for(int i=1; i<=n; ++i) {
            int p = sz-1;
            while(p&&cross(stk[p]-stk[p-1], Pos(1,2*sum[i]))<=0)
                --p;
            dp[i] = stk[p].y+m+sum[i]*(sum[i]-2*stk[p].x);
            Pos now = Pos(sum[i], dp[i]+sum[i]*sum[i]);
            while(sz>1&&cross(stk[sz-1]-stk[sz-2], now-stk[sz-2])<=0)
                --sz;
            stk[sz++] = now;
        }
        printf("%I64d\n", dp[n]);
    }
    return 0;
}
```

* * *

## 插头DP
```c++
//hdu-2167
#include <bits/stdc++.h>//don't use this in poj, fzu, zoj
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = 1e5+1e3;
int n, a[16][16];
int cnt[2][1<<16];
void update(int& a, int b) {
    if(a<b) a=b;
}
bool test(int a, int b) {
    return a&(1<<b);
}
bool valid(int x, int y, int k) {
    if(x==0) return y?!test(k,0):true;
    else {
        if(y==0) return !test(k,n-1)&&!test(k,n-2);
        else if(y==n-1) return !test(k,0)&&!test(k,n)&&!test(k, n-1);
        else return !test(k,0)&&!test(k,n)&&!test(k,n-1)&&!test(k,n-2);
    }
}
void solve() {
    int maxst=1<<(n+1);
    int now=0, pre=1;
    memset(cnt[now], 0, sizeof(int)*maxst);
    for(int i=0 ; i<n; ++i) {
        for(int j=0 ; j<n; ++j) {
            swap(now, pre);
            memset(cnt[now], 0, sizeof(int)*maxst);
            for(int k=0 ; k<maxst; ++k) {
                update(cnt[now][(k<<1)&~maxst], cnt[pre][k]);
                if(valid(i, j, k))
                    update(cnt[now][(k<<1|1)&~maxst], cnt[pre][k]+a[i][j]);
            }
        }
    }
    int ans=0;
    for(int i=0 ; i<maxst; ++i) {
        update(ans, cnt[now][i]);
    }
    pf("%d\n", ans);
}
int input() {
    static char line[MAXN];
    int& curx=n, cury;
    curx=0;
    CLR(a);
    while(gets(line)) {
        cury=0;
        if(line[0]==' '||line[0]==0) return curx;
        int l=strlen(line);
        for(int i=0 ; i<l; ++i) {
            if(line[i]==' ') ++cury;
            else
            (a[curx][cury]*=10)+=line[i]-'0';
        }
        ++curx;
    }
    return curx;
}
void output() {
    for(int i=0 ; i<n; ++i) {
        for(int j=0; j<n; ++j)
            printf("%d ", a[i][j]);
        puts("");
    }
}
int main() {
    int t=1, ca=0;
    while(input()) {
        //output();
        solve();
    }
    return 0;
}
```
***

## 数位DP

```c++
//CF 55D
#include<iostream>
#include<cstring>
using namespace std;
typedef long long ll;
const int MAXL = 20;
const int MAXN = 2530;
const int MAXS = 50;
ll dp[MAXL][MAXN][MAXS];
int d[MAXL];
int hashv[MAXN];
void init() {
    int s=2520, cnt=0;
    for(int i=1; i<50; ++i)
        if(s%i==0) {
            hashv[i] = ++cnt;
            hashv[s/i] = 49-cnt;
        }
}
ll gcd(int a, int b) {
    return b?gcd(b,a%b):a;
}
ll dfs(int pos, int r, int lcm, bool lim) {
    if(pos==-1) {
        if(r%lcm==0) return 1;
        return 0;
    }
    int hv = hashv[lcm];
    if(!lim&&dp[pos][r][hv]!=-1)
        return dp[pos][r][hv];
    int up=lim?d[pos]:9;
    ll ret=0;
    for(int i=0; i<=up; ++i) {
        int nr, nlcm;
        if(i) {
            nr = (r*10+i)%2520;
            nlcm = lcm*i/gcd(lcm,i);
        }
        else {
            nr = r*10%2520;
            nlcm = lcm;
        }
        ret += dfs(pos-1, nr, nlcm, lim&&i==up);
    }
    if(!lim) dp[pos][r][hv]=ret;
    return ret;
}
ll getsum(ll x) {
    int pos=0;
    while(x) {
        d[pos++] = x%10;
        x/=10;
    }
    return dfs(pos-1, 0, 1, true);
}
int main() {
    int t;cin >> t;
    init();
    memset(dp, -1, sizeof(dp));
    while(t--) {
        ll l, r;
        cin >> l >> r;
        cout << getsum(r)-getsum(l-1) <<'\n';
    }
    return 0;
}
```
***

## 大背包
```c++
#include<cstdio>
#include<algorithm>
using namespace std;
typedef long long ll;
typedef pair<ll,ll> Pair;
#define t first
#define w second
const int MAXN = 101;
Pair a[MAXN], sum[MAXN];
int n, m;
ll ans;
void dfs(int tot, ll m, ll get) {
	if(sum[tot].t<=m) {
		ans = max(ans, get+sum[tot].w);
		return;
	} else ans = max(ans, get);
	if(a[1].t>m) return;
	if(get+sum[tot].w<=ans) return;
	dfs(tot-1, m, get);
	if(a[tot].t<=m) dfs(tot-1,m-a[tot].t,a[tot].w+get);
}
int main() {
	while(~scanf("%d%d", &n, &m)) {
		for(int i=1; i<=n; ++i)
			scanf("%I64d%I64d", &a[i].t, &a[i].w);
		sort(a+1, a+n+1);
		sum[0] = make_pair(0ll,0ll);
		for(int i=1; i<=n; ++i)
			sum[i] = make_pair(sum[i-1].t+a[i].t, sum[i-1].w+a[i].w);
		while(a[n].t>m) --n;
		ans = 0;
		dfs(n, m, 0);
		printf("%I64d\n", ans);
	}
	return 0;
}
```

* * *

# 计算几何
## 凸包
```c++
#include<cstdio>
#include<cmath>
#include<algorithm>
using namespace std;
const double Pi = acos(-1.0);
const int MAXN = 1024;
struct Pos{
    int x, y;
    Pos(int x=0, int y=0)
        :x(x), y(y){}
    Pos operator - (const Pos& b)
    {return Pos(x-b.x, y-b.y);}
    bool operator <(const Pos& b)
        const{
            return x==b.x?y<b.y:x<b.x;
    }
    double mod() {
        return sqrt(x*x+y*y+0.0);
    }
}P[MAXN], stk[MAXN];
int n, l, sz;
int Det(const Pos& a, const Pos& b) {
    return a.x*b.y-a.y*b.x;
}
int main() {
    scanf("%d", &n);
    for(int i=0; i<n; ++i)
    scanf("%d%d", &P[i].x, &P[i].y);
    sort(P, P+n);
    sz=0;
    for(int i=0; i<n; ++i) {
        while(sz>1&&Det(stk[sz-1]-stk[sz-2],P[i]-stk[sz-2])<=0)
            sz--;
        stk[sz++] = P[i];
    }
    int k=sz;
    for(int i=n-2; i>-1; --i) {
        while(sz>k&&Det(stk[sz-1]-stk[sz-2],P[i]-stk[sz-2])<=0)
            sz--;
        stk[sz++] = P[i];
    }
    return 0;
}
```

* * *

## 动态凸包
```c++
#include<cstdio>
#include<map>
using namespace std;
typedef long long ll;
#define fi first
#define se second
struct Point{
    int x, y;
    Point(int xx=0, int yy=0)
        :x(xx), y(yy){}
    Point operator - (const Point& b) const{
        return Point(x-b.x, y-b.y);
    }
};
ll cross(const Point& a, const Point& b) {
    return 1ll*a.x*b.y-1ll*a.y*b.x;
}
struct Convex{
    map<int,int> mp[2];
    map<int,int>::iterator l, l2, r, r2;
    Convex() {
        mp[0].clear();//upper
        mp[1].clear();//under
    }
    bool isin(map<int,int>&mp, int x, int y) {
        if(mp.find(x)!=mp.end()) return y>=mp[x];
        l = r = mp.upper_bound(x);
        if(r==mp.begin()||r==mp.end()) return false;
        --l;
        return cross(Point(x-l->fi,y-l->se), Point(r->fi-x,r->se-y))<=0;
    }
    bool isin(int x, int y) {
        return isin(mp[0],x,-y)&&isin(mp[1],x,y);
    }
    void insert(map<int,int>& mp, int x, int y) {
        mp[x] = y;
        r = mp.upper_bound(x);
        if(r!=mp.end()) {
            r2 = r;
            ++r2;
            while(r2!=mp.end()) {
                if(cross(Point(r->fi-x,r->se-y), Point(r2->fi-x, r2->se-y))<=0)
                {
                    mp.erase(r);
                    r = r2;
                }
                else break;
                ++r2;
            }
        }
        l = mp.lower_bound(x);
        if(l==mp.begin()) return;
        --l;
        if(l!=mp.begin()) {
            l2 = l;
            while(true) {
                --l2;
                if(cross(Point(l->fi-l2->fi, l->se-l2->se), Point(x-l->fi, y-l->se))<=0)
                {
                    mp.erase(l);
                    l = l2;
                }
                else break;
                if(l2==mp.begin()) break;
            }
        }
    }
    void insert(int x, int y) {
        if(!isin(mp[0], x, -y))
            insert(mp[0], x, -y);
        if(!isin(mp[1], x, y))
            insert(mp[1], x, y);
    }
    void print() {
        for(int i=0; i<2; ++i) {
            for(l=mp[i].begin(); l!=mp[i].end(); ++l)
                printf("(%d,%d) ", l->fi, l->se*(i?1:-1));
            puts("");
        }
        puts("");
    }
}conv;
```

* * *

## 旋转卡壳法
```c++
int rotating_calipers(Point *ch,int n)
{
    int q=1,ans=0;
    ch[n]=ch[0];
    for(int p=0;p<n;p++)
    {
        while(cross(ch[p+1],ch[q+1],ch[p])>cross(ch[p+1],ch[q],ch[p]))
            q=(q+1)%n;
        ans=max(ans,max(dist2(ch[p],ch[q]),dist2(ch[p+1],ch[q+1])));            
    }
    return ans;
}
```
## 矩形面积求并
```c++
#include <iostream>
#include <algorithm>
 using namespace std;
 const int maxn=110;

 struct LINE{
   double  x, y_down, y_up;
   int  flag;
   bool operator<(const LINE &a)const
   {
      return  x<a.x;
   }
 }line[2*maxn];

 struct TREE{
    double  y_down, y_up;
    double  x;
    int     cover; //用以表示加进线段树中的线段次数
    bool    flag; //此标记用来表示是否有超元线段；为了处理方便加上去的
 }tree[1000*maxn];
 int     n;
 double  x1, y1, x2, y2;
 int     index=0;
 double  y[2*maxn];

 void build(int i, int l, int r) {
        tree[i].x = -1; //-1表示该区间已经没有线段
        tree[i].cover = 0; //表示该区间上有多少条线段；左边线段加进去则++，右边线段加进去则--
        tree[i].y_down = y[l];
        tree[i].y_up = y[r];
        tree[i].flag = false;
        if(l+1==r) {
            tree[i].flag = true; //flag==true表示达到了叶子节点
            return;
        }
        int mid=(l+r)>>1;
        build(2*i, l, mid);
        build(2*i+1, mid, r);
 }
 double insert(int i, double x, double l, double r, int flag) {//flag表示为左边还是右边
     if (r<=tree[i].y_down || l>=tree[i].y_up) return 0;
     if (tree[i].flag)  {
         if (tree[i].cover > 0) {//递归到了叶子节点
              double temp_x = tree[i].x;
              double ans=(x-temp_x)*(tree[i].y_up - tree[i].y_down);
              tree[i].x = x;   //定位上一次的x
              tree[i].cover += flag;
              return ans;
         }
         else {
             tree[i].cover += flag;
             tree[i].x = x;
             return 0;
         }
     }
    double ans1, ans2;
    ans1 = insert(2*i, x, l, r, flag);
    ans2 = insert(2*i+1, x, l, r, flag);
    return ans1+ans2;
 }
 int main( ) {
     int  count=0;
     while (scanf("%d", &n)!=EOF&&n) {
         index = 1;
         for (int i=1; i<=n; i++) {
             scanf("%lf%lf%lf%lf", &x1, &y1, &x2, &y2);
             y[index] = y1;
             line[index].x = x1;
             line[index].y_down = y1;
             line[index].y_up = y2;
             line[index].flag = 1; //1表示左边
             index++;
             y[index] = y2;
             line[index].x = x2;
             line[index].y_down = y1;
             line[index].y_up = y2;
             line[index].flag = -1; //-1表示右边
             index++;
         }
         sort(&y[1], &y[index]); //把所有的纵坐标按从小到大排序,把1写成了0，WA一次
         sort(&line[1], &line[index]);
         build(1, 1, index-1);
         double ans=0;
         for (int i=1;i<index; i++)
             ans += insert(1, line[i].x, line[i].y_down, line[i].y_up, line[i].flag);
         printf("Test case #%d\nTotal explored area: %.2f\n\n", ++count, ans);
     }
     return 0;
}
```

***

## Geo模板

```c++
#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
using namespace std;
const double eps=1e-8;
const double PI=acos(-1.0);
struct Point
{
    double x,y;
    Point() {}
    Point(double a,double b)
    {
        x=a,y=b;
    }
    Point operator + (const Point &b)const
    {
        return Point(x+b.x,y+b.y);
    }
    Point operator - (const Point &b)const
    {
        return Point(x-b.x,y-b.y);
    }
    double operator ^(const Point &b)const
    {
        return x*b.y-y*b.x;
    }
    double operator *(const Point &b)const
    {
        return x*b.x+y*b.y;
    }
    bool operator == (const Point &b)const
    {
        return x==b.x&&y==b.y;
    }
};
struct Line
{
    Point s,e;
    double ang;
    Line(){}
    Line(Point a,Point b)
    {
        s=a,e=b;
        ang=atan2(b.y-a.y,b.x-a.x);
    }
    Point operator &(const Line &b)const
    {
        Point res=s;
        double t=((s-b.s)^(b.s-b.e))/((s-e)^(b.s-b.e));
        res.x+=(e.x-s.x)*t;
        res.y+=(e.y-s.y)*t;
        return res;
    }
};
double get_dis(Point a,Point b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
bool Point_on_segment(Point q,Point s,Point e)
{
    if(fabs((q-s)^(e-s))<eps)
        if(min(s.x,e.x)<=q.x&&q.x<=max(s.x,e.x)&&min(s.y,e.y)<=q.y&&q.y<=max(s.y,e.y))
            return 1;
    return 0;
}
bool Segment_cross(Line a,Line b)
{
    if(max(a.e.x,a.s.x)>=min(b.s.x,b.e.x)&&max(b.s.x,b.e.x)>=min(a.s.x,a.e.x)&&max(a.s.y,a.e.y)>=min(b.s.y,b.e.y)&&max(b.s.y,b.e.y)>=min(a.s.y,a.e.y))
        if(((a.s-b.s)^(b.e-b.s))*((b.e-b.s)^(a.e-b.s))>-eps&&((b.s-a.s)^(a.e-a.s))*((a.e-a.s)^(b.e-a.s))>-eps)
            return 1;
    return 0;
}
bool Segment_Line_cross(Line a,Line b)
{
    return ((a.s-b.s)^(b.e-b.s))*((b.e-b.s)^(a.e-b.s))>-eps;
}
Point Point_to_Line(Point P,Line L)
{
    Point result;
    double t = ((P-L.s)*(L.e-L.s))/((L.e-L.s)*(L.e-L.s));
    result.x = L.s.x + (L.e.x-L.s.x)*t;
    result.y = L.s.y + (L.e.y-L.s.y)*t;
    return result;
}
Point Point_to_segment(Point P,Line L)
{
    Point result;
    double t = ((P-L.s)*(L.e-L.s))/((L.e-L.s)*(L.e-L.s));
    if(t>-eps&&t-1<-eps)
    {
        result.x = L.s.x + (L.e.x-L.s.x)*t;
        result.y = L.s.y + (L.e.y-L.s.y)*t;
    }
    else
        result=get_dis(P,L.s)-get_dis(P,L.e)<-eps?L.s:L.e;
    return result;
}
double Cal_area(Point p[],int n)
{
    double res=0;
    for (int i=1;i<n-1;i++)
        res+=(p[i]-p[0])^(p[i+1]-p[0]);
    return res/2.0;
}
Point Rotate(Point vec,double rad)//Anti-clockwise
{
    return Point(vec.x*cos(rad)-vec.y*sin(rad),vec.x*sin(rad)+vec.y*cos(rad));
}
Point Normal_vector(Point vec)
{
    double l=sqrt(vec*vec);
    /*if(l>-eps&&l<eps)
        puts("Zero vector");*/
    return Point(-vec.y/l,vec.x/l);
}
int Point_in_poly(Point p,Point poly[],int n)
{
    int cnt=0;
    Line ray,side;
    ray.s=p,ray.e.y=p.y,ray.e.x=-INF;
    for (int i=0;i<n;i++)
    {
        side.s=poly[i],side.e=poly[(i+1)%n];
        if(Point_on_segment(p,side.s,side.e))
            return 0;
        if(fabs(side.s.y-side.e.y)<eps)
            continue;
        if(Point_on_segment(side.s,ray.s,ray.e)&&side.s.y-side.e.y>eps)
            cnt++;
        else if(Point_on_segment(side.e,ray.s,ray.e)&&side.e.y-side.s.y>eps)
            cnt++;
        else if(Segment_cross(ray,side))
            cnt++;
    }
    return cnt&1?1:-1;
}
int Segment_in_poly(Line a,Point poly[],int n)
{
    if(Point_in_poly(a.s,poly,n)==0||Point_in_poly(a.e,poly,n))
        return 0;
    for (int i=0;i<n;i++)
    {
        if(Segment_Line_cross(a,Line(poly[i],poly[(i+1)%n])))
            return 0;
        if(Point_on_segment(poly[i],a.s,a.e)&&Segment_Line_cross(Line(poly[(i+1)%n],poly[(i+n-1)%n]),a))
            return 0;
    }
    return 1;
}
bool Is_convex(Point poly[],int n)
{
    bool s[3]={0,0,0};
    for (int i=0;i<n;i++)
    {
        int pos;
        if(fabs((poly[(i+1)%n]-poly[i])^(poly[(i+2)%n]-poly[i]))<eps)
            pos=0;
        else
            pos=((poly[(i+1)%n]-poly[i])^(poly[(i+2)%n]-poly[i]))>eps?1:-1;
        s[pos+1]=1;
        if(s[0]&&s[2])
            return 0;
    }
    return 1;
}
bool Point_cmp(Point a,Point b)
{
    if(a.x==b.x)
        return a.y<b.y;
    return a.x<b.x;
}
int Convex_hull(Point p[],int n,Point ans[])
{
    sort(p,p+n,Point_cmp);
    n=unique(p,p+n)-p;
    int m=0;
    for (int i=0;i<n;i++)
    {
        while(m>1&&((ans[m-1]-ans[m-2])^(p[i]-ans[m-2]))<=0)
            m--;
        ans[m++]=p[i];
    }
    int k=m;
    for (int i=n-2;i>=0;i--)
    {
        while(m>k&&((ans[m-1]-ans[m-2])^(p[i]-ans[m-2]))<=0)
            m--;
        ans[m++]=p[i];
    }
    if(n>1)
        m--;
    return m;
}
double Rotating_Calipers(Point poly[],int n)//Anti-clockwise
{
    int j=1;
    double ans=0;
    poly[n]=poly[0];
    for (int i=0;i<n;i++)
    {
        while(((poly[i]-poly[i+1])^(poly[i]-poly[j]))-((poly[i]-poly[i+1])^(poly[i]-poly[j+1]))<-eps)
            j=(j+1)%n;
        ans=max(ans,max(get_dis(poly[i],poly[j]),get_dis(poly[i+1],poly[j])));
    }
    return ans;
}
bool HPI_cmp(Line a,Line b)
{
    if(fabs(a.ang-b.ang)>eps)
        return a.ang<b.ang;
    return ((a.s-b.s)^(b.e-b.s))<-eps;//������ʱ�뷽���Ǹ�
}
int HPI(Line line[],int n,Point ans[])
{
    Line Q[MAXN];
    int tot=1,cnt=0;
    sort(line,line+n,HPI_cmp);
    for (int i=1;i<n;i++)
        if(fabs(line[i].ang-line[i-1].ang)>eps)
            line[tot++]=line[i];
    int head=0,tail=1;
    Q[0]=line[0];
    Q[1]=line[1];
    for(int i = 2; i < tot; i++)
    {
        /*if(fabs((Q[tail].e-Q[tail].s)^(Q[tail-1].e-Q[tail-1].s))<eps||fabs((Q[head].e-Q[head].s)^(Q[head+1].e-Q[head+1].s))<eps)
            return 0;*/
        while(head < tail && (((Q[tail]&Q[tail-1])-line[i].s)^(line[i].e-line[i].s)) > eps)
            tail--;
        while(head < tail && (((Q[head]&Q[head+1])-line[i].s)^(line[i].e-line[i].s)) > eps)
            head++;
        Q[++tail] = line[i];
    }
    while(head < tail && (((Q[tail]&Q[tail-1])-Q[head].s)^(Q[head].e-Q[head].s)) > eps)
        tail--;
    while(head < tail && (((Q[head]&Q[head-1])-Q[tail].s)^(Q[tail].e-Q[tail].e)) > eps)
        head++;
    if(tail<=head+1)
        return 0;
    for (int i=head;i<tail;i++)
        ans[cnt++]=Q[i]&Q[i+1];
    if(head<tail-1)
        ans[cnt++]=Q[head]&Q[tail];
    return cnt;
}

```

* * *

# 数论
## 离散对数
```c++
int log_mod(int a, int b, int n) {
    int s, v, e;
    s = (int)sqrt(n+0.5);
    v = Inv(pow(a, m, n), n);
    map <int, int> r;
    int cur = 1;
    for(int i=0; i<s; ++i) {
        r[cur] = i;
        cur = cur*a%n;
    }
    for(int i=0; i<s; ++i) {
        if(r.count(b)) return i*s+r[b];
        b = b*v%n;
    }
    return -1;
}
```
## 欧拉函数
```c++
int phi(int n) {
    int m = (int)sqrt(n+0.5);
    int ans = n;
    for(int i=2; i<=m; ++i) if(n%i==0) {
        ans = ans/i*(i-1);
        while(n%i==0) n/=i;
    }
    if(n>1) ans = ans/n*(n-1);
    return ans;
}
```
## 素数测试
```c++
#include<cstdio>
#include<cstdlib>
#include<ctime>
typedef long long ll;
const int MAXN = 5e6;
int a[4] = {2,3,5}, sz;
ll prime[MAXN];
ll pow_mod(ll a, ll x, ll n) {
    ll ret = 1;
    while(x) {
        if(x&1) ret = ret*a%n;
        a = a*a%n;
        x >>=1;
    }
    return ret;
}
bool witness(int a, ll n) {
    ll u=n-1;
    int t = 0;
    while((u&1)==0) {
        u>>=1;
        ++t;
    }
    ll x = pow_mod(a, u, n), px;
    for(int i=0; i<=t; ++i) {
        px = x;
        x = x*x%n;
        if(x==1&&px!=1&&px!=n-1)
            return true;
    }
    if(x!=1) return true;
    return false;
}
bool primetest(ll n) {
    if(n==2||n==3||n==5)
        return true;
    if(n%2==0||n%3==0||n%5==0)
        return false;
    for(int i=0; i<3; ++i)
        if(witness(a[i], n))
            return false;
    return true;
}
```
## 素数个数
```c++
#include<cstdio>
#include<cstring>
#include<cmath>
typedef long long ll;
const int MAXN = 100;
const int MAXM = 10010;
const int MAXP = 777777;
const int MAX = 1000100;
ll dp[MAXN][MAXM];
bool vis[MAX];
int len = 0, prime[MAXP], cnt[MAX];
ll phi(ll m, int n) {
    if(n==0) return m;
    if(prime[n-1]>=m) return 1;
    if(m<MAXM&&n<MAXN) return dp[n][m];
    return phi(m,n-1)-phi(m/prime[n-1],n-1);
}
ll remhel(ll m) {
    if(m<MAX) return cnt[m];
    ll res=0;
    int s, a, c, y;
    s = sqrt(0.9+m);
    y = c = cbrt(0.9+m);
    a = cnt[y];
    res = phi(m,a)+a-1;
    for(int i=a; prime[i]<=s; ++i)
        res = res-remhel(m/prime[i])+remhel(prime[i])-1;
    return res;
}
void init() {
    vis[0] = vis[1] = true;
    for(int i=3; i*i<MAX; i+=2)
        if(!vis[i]) {
            int k=i<<1;
            for(int j=i*i; j<MAX; j+=k)
                vis[j] = true;
        }
    for(int i=1; i<MAX; ++i) {
        cnt[i] = cnt[i-1];
        if((i==2||(i&1))&&!vis[i]) {
            prime[len++] = i;
            ++cnt[i];
        }
    }
    for(int n=0; n<MAXN; ++n)
        for(int m=0; m<MAXM; ++m)
            if(!n) dp[n][m] = m;
            else dp[n][m] = dp[n-1][m]-dp[n-1][m/prime[n-1]];
}
```

## Lucas定理
```c++
ll C(ll x, ll y) {
    if(x>=MOD)
        return C(x/MOD, y/MOD)*C(x%MOD, y%MOD)%MOD;
    else if(x<y) return 0;
    else if(x==y||y==0) return 1;
    else return fac[x]*Inv(fac[y]*fac[x-y]%MOD)%MOD;
}
```

---

## 平方剩余

```c++
#include <bits/stdc++.h>//don't use this in poj, fzu, zoj
using namespace std;
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
#define dbg(x) std::cout << #x"=" << x <<'\n';;
typedef long long ll;
const int MAXN = 1e5+1e3;
struct Node{
    static ll p, omega;
    ll a, b;
    Node(ll aa, ll bb) {
        a=aa%p, b=bb%p;
    }
    Node operator *(const Node& rhs) {
        return Node(a*rhs.a%p+b*rhs.b%p*omega%p, a*rhs.b%p+b*rhs.a%p);
    }
};
ll Node::p, Node:: omega;
ll pow_mod(ll a, ll x, ll p) {
    ll ret=1;
    while (x) {
        if(x&1) ret=ret*a%p;
        a=a*a%p;
        x>>=1;
    }
    if(ret<0) ret+=p;
    return ret;
}
Node pow_mod(Node a, ll x) {
    Node ret(1,0);
    while(x) {
        if(x&1) ret = ret*a;
        a=a*a;
        x>>=1;
    }
    return ret;
}
int lgd(ll a, ll p) {
    if(pow_mod(a, (p-1)>>1, p)==1)
        return 1;
    return -1;
}
ll modsqrt(ll a,ll p){
    a%=p;
    if(p==2||a==0) return a%p;
    if(lgd(a, p)==-1) return -1;
    ll x;
    //if(p%4==3) x=pow_mod(a,(p+1)>>2,p);
    //else{
        ll a_0;
        for(a_0=1; a_0<p; ++a_0) {
            if(lgd(a_0*a_0-a, p)==-1)
                break;
        }
        Node::p=p;
        Node::omega=(a_0*a_0-a)%p;
        Node tmp=Node(a_0, 1);
        tmp = pow_mod(tmp, (p+1)>>1);
        x=tmp.a;
        dbg(tmp.b);
    //}
    if(x<0) x+=p;
    if(x*2>p) x=p-x;
    return x;
}
```
---

# 区间问题
## 莫队算法
```c++
#include<cstdio>
#include<cmath>
#include<algorithm>
#include<cstring>
const int MAXN = 2e5+1e4;
const int MAXV = 1e6+1e3;
typedef long long ll;
using std::sort;
struct Ask{
    int l, r, id, b;
}q[MAXN];
```

---

```c++
int n, t, a[MAXN], sqn;
int cnt[MAXV];
ll ans[MAXN], nv;
bool cmp(const Ask& a, const Ask& b) {
    return a.b==b.b?a.r<b.r:a.b<b.b;
}
void remove(int pos) {
    nv -= (2*cnt[a[pos]]-1ll)*a[pos];
    cnt[a[pos]]--;
}
void add(int pos) {
    nv += (2*cnt[a[pos]]+1ll)*a[pos];
    cnt[a[pos]]++;
}
int main() {
    scanf("%d%d", &n, &t);
    for(int i=1; i<=n; ++i)
        scanf("%d", a+i);
    sqn = (int)sqrt(n+0.5);
    for(int i=0; i<t; ++i) {
        scanf("%d%d", &q[i].l, &q[i].r);
        q[i].id = i;
        q[i].b = q[i].l/sqn;
    }
    sort(q, q+t, cmp);
    int nl=0, nr=0;
    nv = 0;
    for(int i=0; i<t; ++i) {
        while(nr<q[i].r)
            add(++nr);
        while(nr>q[i].r)
            remove(nr--);
        while(nl<q[i].l)
            remove(nl++);
        while(nl>q[i].l)
            add(--nl);
        ans[q[i].id] = nv;
    }
    for(int i=0; i<t; ++i)
        printf("%I64d\n", ans[i]);
    return 0;
}
```
## fzu2224
```c++
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
const int MAXN = 1e4+1e3;
const int MAXK = 20;
typedef long long ll;
struct Que{
    int l, r , id;
    bool operator <(const Que& b) const{
        return r<b.r;
    }
}q[MAXN*10];
struct Seg{
    int l, r, v;
    bool operator <(const Seg& b) const{
        return r<b.r;
    }
}seg[MAXN*MAXK];
```

---

```c++
int n, m, tot, mat[MAXN*MAXK];
int d[MAXN][MAXK], hsh[MAXN*MAXK];
int ans[MAXN*10], tree[MAXN];
int gcd(int a, int b) {
    return b?gcd(b,a%b):a;
}
void add(int pos, int v) {
    for(int i=pos; i<=n; i+=i&(-i))
        tree[i] += v;
}
ll sum(int pos) {
    ll ret = 0;
    for(int i=pos; i>0; i-=i&(-i))
        ret += tree[i];
    return ret;
}
int erfen(int pos, int v) {
    int l=1, r=pos;
    while(l<r) {
        int mid=(l+r)/2;
        int len = pos-mid+1;
        int k=0, cur=0;
        while(1<<(k+1)<=len) ++k;
        cur = gcd(d[pos][k], d[mid+(1<<k)-1][k]);
        if(cur<v) l= mid+1;
        else r = mid;
    }
    return l;
}
int main() {
int T;scanf("%d", &T);
    while(T--)
    {
scanf("%d%d", &n, &m);
        for(int i=1; i<n+1; ++i)
            scanf("%d", d[i]);
        for(int k=1; (1<<k)<n+1; ++k)
            for(int i=n; (1<<k)<i+1; --i)
                d[i][k] = gcd(d[i][k-1],d[i-(1<<k-1)][k-1]);
        for(int i=1; i<m+1; ++i) {
            scanf("%d%d", &q[i].l, &q[i].r);
            q[i].id = i;
        }
        int cnt = 0;
        tot = 0;
        for(int i=1; i<=n; ++i) {
            int now = 0;
            for(int j=i; j>0; --j) {
                now = gcd(now, d[j][0]);
                ++cnt;
                seg[cnt] = (Seg){j,i,now};
                hsh[++tot] = now;
                j = erfen(i,now);
            }
        }
        sort(hsh+1, hsh+1+tot);
        tot=unique(hsh+1,hsh+1+tot)-hsh-1;
        sort(q+1, q+m+1);
        sort(seg+1, seg+1+cnt);
        memset(mat, 0, sizeof(mat));
        memset(tree, 0, sizeof(tree));
        for(int i=1,now=1; i<=m; ++i) {
            for(;now<=cnt&&seg[now].r<=q[i].r; ++now) {
                int pos = lower_bound(hsh+1, hsh+1+tot, seg[now].v)-hsh;
                if(seg[now].l>mat[pos]) {
                    if(mat[pos]) add(mat[pos],-1);
                    add(seg[now].l, 1);
                    mat[pos]=seg[now].l;
                }
            }
            ans[q[i].id] = sum(q[i].r)-sum(q[i].l-1);
        }
        for(int i=1; i<=m; ++i)
            printf("%d\n", ans[i]);
    }
    return 0;
}
```
## Treap
```c++
#include<bits/stdc++.h>
using namespace std;
struct Node{
    Node *ch[2]; //l:0, r:1;
    int r, v, s;// 堆值，小根堆 节点值, 结点数
    Node(int v): v(v) {r = rand(); s=1; ch[0]=ch[1]=NULL;}
    int cmp(int x) const{
        if(x==v) return -1;
        return x<v?0:1;
    }
    void maintain() {
        s = 1;
        if(ch[0]!=NULL) s+= ch[0]->s;
        if(ch[1]!=NULL) s+= ch[1]->s;
    }
};
```

---

```c++
struct Treap {
    Node *root;
    void rotate(Node*& o, int d)//d=0 to left, d=1 to right
    {
        Node* k = o->ch[1^d];
        o->ch[1^d] = k->ch[d];
        k->ch[d] = o;
        o->maintain();
        k->maintain();
        o = k;
    }
    bool find(Node* o, int x) {
        while(o!=NULL) {
            int d = o->cmp(x);
            if(d==-1) return true;
            else o = o->ch[d];
        }
        return false;
    }
    void insert(Node* &o, int x) {
        if(o==NULL)
            o = new Node(x);
        else {
            int d = o->cmp(x);
            if(d!=-1) {
                insert(o->ch[d], x);
                if(o->ch[d]->r > o->r) rotate(o, 1^d);
            }
        }
    }
    void remove(Node* &o, int x) {
        int d = o->cmp(x);
        if(d==-1) {
            if(o->ch[0]==NULL)
                o = o->ch[1];
            else if(o->ch[1]==NULL)
                o = o->ch[0];   
            else {
                int d2 = (o->ch[0]->r > o->ch[1]->r)?1:0;
                rotate(o, d2);
                remove(o->ch[d2], x);
            }
        }
        else remove(o->ch[d], x);
    }
    void remove(int x) {
        if(find(root, x)) remove(root, x);
    }
    void clear(Node*& o) {
        if(o->ch[0]!=NULL) clear(o->ch[0]);
        if(o->ch[1]!=NULL) clear(o->ch[1]);
        delete o;
        o = NULL;
    }
    int Kth(Node* o, int k) {//第k小
        if(o==NULL|| k<=0 ||k > o->s)
            return 0;
        int s = (o->ch[0]==NULL?0:o->ch[0]->s);
        if(k==s+1) return o->v;
        if(k<s) return Kth(o->ch[0], k);
        else return Kth(o->ch[1], k-s-1);
    }
};
```
## 树上莫队cot2
```c++
#include<cstdio>
#include<cmath>
#include<cstring>
#include<algorithm>
using namespace std;
const int MAXN = 41000;
const int MAXM = 1e5+1e3;
const int MAXL = 18;
int anc[MAXN][MAXL], d[MAXN], dfn[MAXN];
int n, m, w[MAXN], sqn, ne, dfs_clock;
int ans[MAXM], b[MAXN], cnt[MAXN];
bool vis[MAXN];
struct Edge{
	int to;
	Edge* next;
}e[MAXN*2], *head[MAXN];
struct Que{
	int id, u, v;
	bool operator < (const Que& b) const {
		int p1=dfn[u]/sqn, p2=dfn[b.u]/sqn;
		return p1==p2?dfn[v]<dfn[b.v]:p1<p2;
	}
}q[MAXM];
struct Set{
	int res;
	Set() {
		res = 0;
	}
	void Xor(int pos) {
		if(vis[pos]) {
			vis[pos] = false;
			if(--cnt[w[pos]]==0)
				--res;
		}
		else {
			vis[pos] = true;
			if(cnt[w[pos]]++==0)
				++res;
		}
	}
}cup;
void dfs(int o, int pa) {
	anc[o][0] = pa;
	d[o] = d[pa]+1;
	dfn[o] = ++dfs_clock;
	for(Edge* p=head[o]; p; p=p->next)
		if(p->to!=pa) dfs(p->to, o);
}
void lca_init() {
	dfs(1,1);
	for(int j=1; (1<<j)<=n; ++j)
		for(int i=1; i<=n; ++i)
			anc[i][j] = anc[anc[i][j-1]][j-1];
}
int LCA(int x, int y) {
	if(d[x]<d[y]) swap(x,y);
	int log=1;
	while((1<<log)<=d[x]) ++log;
	--log;
	for(int i=log; i>=0; --i)
		if(d[x]-(1<<i)>=d[y]) x=anc[x][i];
	while(x!=y) {
		int i=0;
		while(i<MAXL&&anc[x][i]!=anc[y][i])
			++i;
		if(i==0) return anc[x][0];
		x = anc[x][i-1];
		y = anc[y][i-1];
	}
	return x;
}
void Modify(int u, int v) {
	int lca = LCA(u,v);
	while(u!=lca) {
		cup.Xor(u);
		u = anc[u][0];
	}
	while(v!=lca) {
		cup.Xor(v);
		v = anc[v][0];
	}
}
void Mo() {
	int l=q[0].u, r=q[0].u;
	sort(q, q+m);
	for(int i=0; i<m; ++i) {
		if(l!=q[i].u) {
			Modify(l, q[i].u);
			l = q[i].u;
		}
		if(r!=q[i].v) {
			Modify(r,q[i].v);
			r = q[i].v;
		}
		int lca = LCA(l, r);
		cup.Xor(lca);
		ans[q[i].id] = cup.res;
		cup.Xor(lca);
	}
}
void addedge(int fr, int to) {
	e[ne].to = to;
	e[ne].next = head[fr];
	head[fr] = e+ne++;
}
int main() {
	scanf("%d%d", &n, &m);
	sqn = (int)sqrt(n+0.5)+1;
	for(int i=1; i<=n; ++i) {
		scanf("%d", w+i);
		b[i] = w[i];
	}
	sort(b+1, b+n+1);
	int c= unique(b+1, b+n+1)-b-1;
	for(int i=1; i<=n; ++i)
		w[i] = lower_bound(b+1, b+c+1, w[i])-b-1;
	int u, v;
	for(int i=1; i<n; ++i) {
		scanf("%d%d", &u, &v);
		addedge(u,v);
		addedge(v,u);
	}
	lca_init();
	for(int i=0; i<m; ++i) {
		q[i].id = i;
		scanf("%d%d", &u, &v);
		if(u>v) swap(u,v);
		q[i].u = u, q[i].v = v;
	}
	Mo();
	for(int i=0; i<m; ++i)
		printf("%d\n", ans[i]);
	return 0;
}
```
## Splay
```c++
#include<bits/stdc++.h>
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define sf scanf
#define pf printf
#define lch(x) ch[x][0]
#define rch(x) ch[x][1]
typedef long long ll;
using namespace std;
int n, m;
const int MAXN = 1e5+1e3;
    int f[MAXN], ch[MAXN][2], key[MAXN];
    int cnt[MAXN], size[MAXN], sz, root;
    int flip[MAXN];
    void clear(int x) {
        lch(x)=rch(x)=f[x]=cnt[x]=key[x]=size[x]=0;
        flip[x]=0;//
    }
    int newnode(int v, int fa) {
        ++sz;
        lch(sz)=rch(sz)=0;
        cnt[sz]=size[sz]=1;
        key[sz]=v, f[sz]=fa;
        flip[sz]=0;//
        return sz;
    }
    void push_down(int o) {
        if(flip[o]) {
            flip[o]=0;
            swap(lch(o), rch(o));
            if(lch(o)) flip[lch(o)]^=1;
            if(rch(o)) flip[rch(o)]^=1;
        }
    }
    void travel(int o) {
        push_down(o);
        if(lch(o)) travel(lch(o));
        printf("%d\n", key[o]);
        if(rch(o)) travel(rch(o));
    }
    int get(int x) {
        return ch[f[x]][1]==x;//0为左儿子，1为右儿子
    }
    void update(int x) {
        if(x) {
            size[x]=cnt[x];
            if(lch(x)) size[x]+=size[lch(x)];//默认x==0代表x是空节点
            if(rch(x)) size[x]+=size[rch(x)];
        }
    }
    void rotate(int x) {//如果他是左儿子则右旋，否则反之
        int old=f[x], oldf=f[old], which=get(x);
        ch[old][which] = ch[x][1^which];
        if(ch[old][which]) f[ch[old][which]] = old;
        ch[x][1^which]=old;
        f[old]=x;
        f[x]=oldf;
        if(oldf) ch[oldf][ch[oldf][1]==old]=x;
        update(old);// update 的顺序不要搞错了
        update(x);
        if(old==root) root=x;
    }
    void splay(int x, int anc) {
        for(int fa; (fa=f[x])!=anc; rotate(x))
            if(f[fa]!=anc)
                rotate(get(fa)==get(x)?fa:x);
    }
    void insert(int v) {
        if(root==0) {
            root = newnode(v,0);
            return;
        }
        int now=root;
        while(key[now]!=v&&ch[now][key[now]<v])
            now = ch[now][key[now]<v];
        if(key[now]==v) {
            ++cnt[now];
            ++size[now];
            splay(now,0);
        }
        else {
            int u= newnode(v,now);
            ch[now][key[now]<v]=u;
            update(now);
            splay(u,0);
        }
    }
    int pre(int o) {
        int now=lch(o);
        while(rch(now)) now=rch(now);
        return now;
    }
    int next(int o) {
        int now=rch(o);
        while(lch(now)) now=lch(now);
        return now;
    }
    void delnode(int o) {
        splay(o, 0);
        if(cnt[root]>1) {
            --cnt[root], --size[root];
            return;
        }
        if(!lch(root)&&!rch(root)) {
            clear(root);
            root=0;
            return;
        }
        if(!lch(root)) {
            root=rch(root);
            clear(f[root]);
            f[root]=0;
        }
        else if(!rch(root)) {
            root=lch(root);
            clear(f[root]);
            f[root]=0;
        }
        else {
            int leftbig=pre(root), oldroot=root;
            splay(leftbig, 0);
            ch[root][1]=ch[oldroot][1];
            f[ch[oldroot][1]] = root;
            clear(oldroot);
            update(root);
        }
    }
    int Kth(int k, int& o) {
        int now=o;
        while(k) {
            push_down(now);
            if(size[lch(now)]>=k)
                now=lch(now);
            else {
                k-=size[lch(now)];
                if(k<=cnt[now]) {
                    k=0;
                    break;
                }
                k-=cnt[now];
                now=rch(now);
            }
        }
        splay(now, f[o]);
        o=now;
        return now;
    }
    void split(int& o, int& r, int k) {
        if(!k) {
            r=o;
            o=0;
            return;
        }
        Kth(k,o);
        r=rch(o);
        f[r]=0;
        rch(o)=0;
    }
    void PT(int a, int b) {
        int l, m, r;
        l=root;
        split(l, m, a-1);
        split(m, r, b-a+1);
        if(l==0) l=r;
        else {
            rch(l)=r;
            f[r]=l;
            update(l);
        }
        rch(m)=lch(m);
        lch(m)=l;
        if(l) f[l]=m;
        if(rch(m)) flip[rch(m)]^=1;
        update(m);
        root = m;
    }
int main() {
    int m;
    sf("%d%d", &n, &m);
    rep(i, 1, n+1)
        insert(i);
    while(m--) {
        int a, b;
        sf("%d%d", &a, &b);
        PT(a,b);
    }
    travel(root);

    return 0;
}

```

***

## HJTree
```c++
// HDU 2665
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
const int MAXN = 1e5+10;
const int MAXM = MAXN*20;
struct Node {
    int ch[2], cnt;
}tr[MAXM];
int a[MAXN],b[MAXN], tot, rt[MAXN];
void update(int p, int& o, int v, int l, int r) {
    tr[o=++tot] = tr[p];
    ++tr[o].cnt;
    if(l < r) {
        int m = (l+r)>>1;
        if(v>m) {
            tr[o].ch[0] = tr[p].ch[0];
            update(tr[p].ch[1], tr[o].ch[1], v, m+1, r);
        } else {
            tr[o].ch[1] = tr[p].ch[1];
            update(tr[p].ch[0], tr[o].ch[0], v, l, m);
        }
    }
}
int ask(int p, int o, int k, int l, int r) {
    if(l==r) return l;
    else {
        int lsz = tr[tr[o].ch[0]].cnt-tr[tr[p].ch[0]].cnt;
        int m = (l+r)>>1;
        if(k>lsz) return ask(tr[p].ch[1], tr[o].ch[1], k-lsz,m+1, r);
        else return ask(tr[p].ch[0], tr[o].ch[0], k, l, m);
    }
}
int main() {
    int T; scanf("%d", &T);
    while(T--) {
        int n, m, sz;
        scanf("%d%d", &n, &m);
        for(int i=1; i<=n; ++i) {
            scanf("%d", a+i);
            b[i] = a[i];
        }
        sort(b+1, b+n+1);
        sz = unique(b+1, b+n+1)-b-1;
        for(int i=1; i<=n; ++i)
            a[i] = lower_bound(b+1, b+sz+1, a[i])-b;
        tot = 0;
        for(int i=1; i<=n; ++i)
            update(rt[i-1], rt[i], a[i], 1, sz);
        while(m--) {
            int l, r, k;
            scanf("%d%d%d", &l, &r, &k);
            printf("%d\n", b[ask(rt[l-1], rt[r], k, 1, sz)]);
        }
    }
    return 0;
}
```
***

## 划分树
```c++
#include<stdio.h>
#include<algorithm>
using namespace std;
#define N 100100
#define rep(a,b,c) for(int a=b; a<c; ++a)
int data[N];

struct node
{
    int v[N];
    int num[N];
} td[31];
void build(int l,int r,int dep)
{
    if(l>=r)return;
    int i,mid=(l+r)>>1,midd=data[mid],ant=mid-l+1;
    // ant保存有多少和sorted[mid]一样大的数进入左孩子
    int ln=l-1,rn=mid;
    for(i=l; i<=r; i++)
        if(td[dep].v[i]<midd)
            ant--;
    for(i=l; i<=r; i++)
    {
        if(i==l)
            td[dep].num[i]=0;
        else
            td[dep].num[i]=td[dep].num[i-1];
        if(td[dep].v[i]<midd||(td[dep].v[i]==midd&&ant))
        {
            if(td[dep].v[i]==midd)--ant;
            td[dep+1].v[++ln]=td[dep].v[i];
            ++td[dep].num[i];
        }
        else td[dep+1].v[++rn]=td[dep].v[i];
    }
    build(l,mid,dep+1);
    build(mid+1,r,dep+1);
}
int query(int a,int b,int k,int l,int r,int dep)
{
    int mid=(l+r)>>1;
    if(a==b)return td[dep].v[a];
    int lx,rx=td[dep].num[b],sub,sl,sr;
    lx=(a<=l)?0:td[dep].num[a-1];
    sub=rx-lx;
    if(sub>=k)
        return query(l+lx,l+rx-1,k,l,mid,dep+1);
    else
    {
        sl=a-l-lx;
        sr=b-l-rx;
        return query(mid+1+sl,mid+1+sr,k-sub,mid+1,r,dep+1);
    }
}
int n;
void dbg() {
    rep(i, 0, 6) {
        printf("%d\n", i);
        rep(j, 1, n+1)
            printf("%d%c", td[i].num[j], j==n?'\n':' ');
    }
    rep(i, 0, 6) {
        printf("%d\n", i);
        rep(j, 1, n+1)
            printf("%d%c", td[i].v[j], j==n?'\n':' ');
    }
}
int main()
{
    int i,cas=1,m,l,r,k;
  //  scanf("%d",&cas);
    while(cas--)
    {
        scanf("%d%d",&n,&m);
        for(i=1; i<=n; i++)
        {
            scanf("%d",&td[1].v[i]);
            data[i]=td[1].v[i];
        }
        sort(data+1,data+n+1);
        build(1,n,1);
        dbg();
        m=0;
        while(m--)
        {
            scanf("%d%d%d",&l,&r,&k);
            printf("%d\n",query(l,r,k,1,n,1));
        }
    }
    return 0;
}
```
***
## TreeinTree
```c++
//HYSBZ - 3295
//树套树的空间一般很大，由于树状数组的结点大小一般较小，可以使用动态线段树，当然也可以使用平衡树
#include<cstdio>
#include<cstring>
#include<bitset>
using namespace std;
const int MAXN = 1e5+1e3;
const int MAXM = 9e6;
struct Node{
    int sz, ch[2];
}tr[MAXM];
struct Act{
    int pos, v;
    long long ans;
}act[MAXN];
int rt[MAXN], n, m, sz;
int a[MAXN], id[MAXN];
bitset<MAXN> isdel;
void Add(int& o, int l, int r, int v) {
    if(!o) o=++sz;
    ++tr[o].sz;
    if(l<r) {
        int m=(l+r)/2;
        if(v<=m) Add(tr[o].ch[0], l, m, v);
        else Add(tr[o].ch[1], m+1, r, v);
    }
}
void Add(int pos, int v) {
    while(pos<=n) {
        Add(rt[pos], 1, n, v);
        pos += pos&(-pos);
    }
}
int Ask(int o, int l, int r, int v) {
    if(!tr[o].sz) return 0;
    if(r<=v) return tr[o].sz;
    if(l>v) return 0;
    int m = (l+r)/2;
    int ret=0;
    if(m<=v) ret = tr[tr[o].ch[0]].sz+Ask(tr[o].ch[1], m+1, r, v);
    else ret = Ask(tr[o].ch[0], l, m, v);
    return ret;
}
int Ask(int pos, int v) {
    int ret=0;
    while(pos) {
        ret += Ask(rt[pos], 1, n, v);
        pos -= pos&(-pos);
    }
    return ret;
}
int Asksz(int pos) {
    int ret=0;
    while(pos) {
        ret += tr[rt[pos]].sz;
        pos -= pos&(-pos);
    }
    return ret;
}
int main() {
    scanf("%d%d", &n, &m);
    for(int i=1; i<=n; ++i) {
        scanf("%d", a+i);
        id[a[i]] = i;
    }
    long long res = 0;
    for(int i=1; i<=m; ++i) {
        int x; scanf("%d", &x);
        act[i] = (Act){id[x], x,0ll};
        isdel[id[x]] = true;
    }
    for(int i=n; i; --i)
        if(!isdel[i]) {
            res += Ask(n, a[i]);
            Add(i, a[i]);
        }
    for(int i=m; i; --i) {
        int pos = act[i].pos, v = act[i].v;
        res += Ask(n, v)+Asksz(pos)-2*Ask(pos, v);
        act[i].ans = res;
        Add(pos, v);
    }
    for(int i=1; i<=m; ++i)
        printf("%lld\n", act[i].ans);
    return 0;
}
```

***

# 数学
## FFT
```c++
#include <bits/stdc++.h>//don't use this in poj, fzu, zoj
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = (1<<18)+10;
const double pi = acos(-1.0);

struct Complex{
    double r, i;
    Complex(double _r=0, double _i=0)
        :r(_r), i(_i){}
    Complex operator + (const Complex& b) {
        return Complex(r+b.r, i+b.i);
    }
    Complex operator - (const Complex& b) {
        return Complex(r-b.r, i-b.i);
    }
    Complex operator * (const Complex& b) {
        return Complex(r*b.r-i*b.i, r*b.i+i*b.r);
    }
    Complex operator / (double b) {
        return Complex(r/b, i/b);
    }
    Complex conj() {
        return Complex(r, -i);
    }
};
char s[MAXN], t[MAXN];
int rev[MAXN];
Complex w[MAXN];
int geth(int x) {
    int ret=0;
    while ((1<<ret) < x)
        ++ret;
    return ret;
}
void fft_prepare(int n, int h) {
    rev[0]=0;
    rep(i, 1, n) {
        rev[i]=rev[i>>1]>>1|((i&1)<<(h-1));
    }
    rep(i, 0, n) {
        w[i]=Complex(cos(2*pi*i/n), sin(2*pi*i/n));
        //w[i+(n>>1)]=Complex(-w[i].r, -w[i].i);
    }
    w[n]=Complex(1,0);
}
void dft(Complex *a, int n, int h, int on) {
    rep(i, 0, n) {
        if(i<rev[i]) swap(a[i], a[rev[i]]);
    }
    for(int s=1 ; s<h+1; ++s) {
        int m=1<<s;
        int mid=m>>1;
        int stp=n/m;
        for(Complex *p=a; p<a+n; p+=m) {
            for(int i=0; i<mid; ++i) {
                Complex t=p[mid+i];
                if (on==1) t=t*w[i*stp];
                else t=t*w[n-i*stp];
                p[mid+i] = p[i]-t;
                p[i] = p[i]+t;
            }
        }
    }
    if (on==-1) {
        rep(i, 0, n) {
            a[i]=a[i]/n;
        }
    }
}
void solve() {
    static Complex a[MAXN], b[MAXN];
    static int out[MAXN];
    int n=strlen(s);
    int m=strlen(t);
    int h=geth(n+m);
    int tn=1<<h;
    fft_prepare(tn,h);
    rep(i, 0, tn) {
        b[i]=Complex();
    }
    rep(i, 0, n) {
        b[i].r = s[n-i-1]-'0';
    }
    rep(i, 0, m) {
        b[i].i = t[m-i-1]-'0';
    }
    dft(b, tn, h, 1);
    rep(i, 0, tn) {
        int j=(tn-i)&(tn-1);
        a[i] = (b[i]*b[i]-b[j].conj()*b[j].conj())*Complex(0, -0.25);
    }
    dft(a, tn, h, -1);
    int len=0, x=0;
    CLR(out);
    rep(i, 0, tn) {
        x=int(x+a[i].r+0.5);
        out[i]=x%10;
        x/=10;
        if(out[i]) len=i+1;
    }
    while (x) {
        out[len++]=x%10;
        x/=10;
    }
    while (len>1&&out[len-1]==0) {
        --len;
    }
    if (!len) len=1;
    drep(i, len-1, -1)
        printf("%d", out[i]);
    puts("");
}
int main() {
    while(~sf("%s%s", s, t)) {
        solve();
    }
    return 0;
}
```

---

## NTT
```c++
#include <bits/stdc++.h>
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = 1<<19;
const ll MOD = 998244353ll;
ll inv[MAXN];
ll fact[MAXN], ifact[MAXN];
void init() {
    inv[1]=1;
    rep(i, 2, MAXN)
        inv[i] = -(MOD/i)*inv[MOD%i]%MOD;
    fact[0] = ifact[0] = 1;
    rep(i, 1, MAXN) {
        fact[i] = i*fact[i-1]%MOD;
        ifact[i] = inv[i]*ifact[i-1]%MOD;
    }
}
ll c[MAXN], s[MAXN];
ll mypow(ll a, ll x) {
    ll ret=1;
    while(x) {
        if(x&1) ret=ret*a%MOD;
        a=a*a%MOD;
        x >>= 1;
    }
    return ret;
}
int rev(int x, int n) {
    int ret=0;
    while(n>1) {
        ret =ret<<1|(x&1);
        n>>=1;
        x>>=1;
    }
    return ret;
}
void ntt(ll* a, int n, int t) {
    rep(i, 0, n) {
        int rv = rev(i,n);
        if(i<rv) swap(a[i], a[rv]);
    }
    ll g=1;
    if(t==1) g=3;
    else g=inv[3];
    for(int m=2; m<n+1; m<<=1) {
        ll wm= mypow(g, (MOD-1)/m);
        int mid=m>>1;
        for(ll *p=a; p<a+n; p+=m) {
            ll w=1;
            ll *a1=p, *a2=p+mid;
            rep(i, 0, mid) {
                ll t=w*(*a2);
                *a2 = (*a1-t)%MOD;
                *a1 = (*a1+t)%MOD;
                w = w*wm%MOD;
                ++a1, ++a2;
            }
        }
    }
}
void fft(ll *a, ll* b, int n) {
    int tn=MAXN;
    while(tn/4>=n) tn>>=1;
    ntt(a, tn, 1);
    ntt(b, tn, 1);
    rep(i, 0, tn)
        a[i] = a[i]*b[i]%MOD;
    ntt(a, tn, -1);
    rep(i, 0, tn)
        a[i]=a[i]*inv[tn]%MOD;
}
int main() {
    int n;
    init();
    while(~sf("%d", &n)) {
        ++n;
        CLR(s); CLR(c);
        drep(i, n-1, -1) {
            sf("%I64d", c+i);
            c[i] = c[i]*fact[n-1-i]%MOD;
        }
        int m;
        sf("%d", &m);
        ll sum=0, tot=1;
        while(m--) {
            ll x;sf("%I64d", &x);
            sum = sum-x;
        }
        sum %= MOD;
        rep(i, 0, n) {
            s[i] = tot*ifact[i]%MOD;
            tot = tot*sum%MOD;
        }
        fft(c,s,n);
        rep(i, 0, n) {
            ll ans=c[n-1-i]*ifact[i]%MOD;
            if(ans<0) ans+=MOD;
            printf("%lld ", ans);
```
---

## FWT
```c++
#include <bits/stdc++.h>//don't use this in poj, fzu, zoj
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = 3e5+1e3;
int a[MAXN], b[MAXN], n;
void fwt(int* a, int n, int t) {
    for(int d=2; d<=n; d<<=1) {
        int m=d>>1;
        for(int i=0; i<n; i+=d) {
            for(int j=0; j<m; ++j) {
                int t=a[i+j+m];
                a[i+j+m]=a[i+j]-t;
                a[i+j] += t;
            }
        }
    }
    if(t==-1) {
        rep(i, 0, n) {
            a[i]/=n;
        }
    }
}
void solve() {
    CLR(a);
    CLR(b);
    sf("%d", &n);
    rep(i, 0, n) {
        sf("%d", a+i);
    }
    rep(i, 0, n) {
        sf("%d", b+i);
    }
    int tn=1;
    while (tn<n)   tn<<=1;
    fwt(a, tn, 1);
    fwt(b, tn, 1);
    rep(i, 0, tn) {
        a[i]*=b[i];
    }
    fwt(a, tn, -1);
    rep(i, 0, n) {
        pf("%d%c", a[i], " \n"[i==n-1]);
    }
}
int main() {
    int t=1, ca=0;
    while(t--) {
        solve();
    }
    return 0;
}

```
---

## 拆系数FFT

```c++
#include <bits/stdc++.h> //don't use this in poj, fzu, zoj
#define rep(a, b, c) for (int(a) = (b); (a) < (c); ++(a))
#define drep(a, b, c) for (int(a) = (b); (a) > (c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = (1 << 17) + 10;
const int MOD = 1e9 + 7;
const double pi = acos(-1.0);
struct Complex {
  double r, i;
  Complex(double _r = 0, double _i = 0) : r(_r), i(_i) {}
  Complex operator+(const Complex &b) { return Complex(r + b.r, i + b.i); }
  Complex operator-(const Complex &b) { return Complex(r - b.r, i - b.i); }
  Complex operator*(const Complex &b) {
    return Complex(r * b.r - i * b.i, r * b.i + i * b.r);
  }
  Complex conj() {
      return Complex(r, -i);
  }
};
ll a[MAXN], b[MAXN], inv[MAXN], c[MAXN];
int geth(int x) {
  int ret = 0, tot = 1;
  while (tot < x) {
    tot <<= 1;
    ++ret;
  }
  return ret;
}
Complex w[MAXN];
int rev[MAXN];
void fft_prepare(int n, int h) {
  rep(i, 0, n >> 1) {
    w[i] = Complex(cos(2 * pi * i / n), sin(2 * pi * i / n));
    w[i + (n >> 1)] = Complex(-w[i].r, -w[i].i);
  }
  w[n] = w[0];
  rev[0] = 0;
  rep(i, 0, n) {
    rev[i] = rev[i >> 1] >> 1;
    if (i & 1)
      rev[i] |= 1 << (h - 1);
  }
}
void dft(Complex *a, int n, int h, int on) {
  rep(i, 0, n) {
    if (i < rev[i])
      swap(a[i], a[rev[i]]);
  }
  for (int m = 2; m < n + 1; m *= 2) {
    int stp = n / m, mid = m >> 1;
    Complex t;
    for (Complex *p = a; p < a + n; p += m) {
      rep(i, 0, mid) {
        if (on > 0)
          t = w[i * stp] * p[mid + i];
        else
          t = w[(n - i * stp) & (n - 1)] * p[mid + i];
        p[mid + i] = p[i] - t;
        p[i] = p[i] + t;
      }
    }
  }
  if (on < 0) {
    rep(i, 0, n) a[i] = Complex(a[i].r / n, a[i].i / n);
  }
}
void conv(ll *a, ll *b, ll* z, int n) {
  int h = geth(n * 2);
  int tn = 1 << h;
  fft_prepare(tn, h);
  static Complex wa[MAXN], wb[MAXN];
  static Complex wc[MAXN], wd[MAXN];
  rep(i, 0, tn) {
    if (i < n) {
      wa[i] = Complex(a[i] >> 15, a[i] & 32767);
      wb[i] = Complex(b[i] >> 15, b[i] & 32767);
    } else
      wb[i] = wa[i] = Complex();
  }
  dft(wa, tn, h, 1);
  dft(wb, tn, h, 1);
  rep(i, 0, tn) {
      int j=(tn-i)&(tn-1);
      Complex dfta = (wa[i]+wa[j].conj())*Complex(0.5,0);
      Complex dftx = (wa[i]-wa[j].conj())*Complex(0, -0.5);
      Complex dftb = (wb[i]+wb[j].conj())*Complex(0.5,0);
      Complex dfty = (wb[i]-wb[j].conj())*Complex(0, -0.5);
      wc[i] = dfta*dftb+dfta*dfty*Complex(0,1);
      wd[i] = dftb*dftx+dftx*dfty*Complex(0,1);
  }
  dft(wc, tn, h, -1);
  dft(wd, tn, h, -1);
  rep(i, 0, tn) {
      ll a=ll(wc[i].r+0.5)%MOD;
      ll b=ll(wc[i].i+0.5)%MOD;
      ll c=ll(wd[i].r+0.5)%MOD;
      ll d=ll(wd[i].i+0.5)%MOD;
      z[i] = ((a<<30)+((b+c)<<15)+d)%MOD;
      if(z[i]<0) z[i] += MOD;
  }
}
void solve() {
  int n, k;
  sf("%d%d", &n, &k);
  rep(i, 0, n) { sf("%I64d", a + i); }
  b[0] = inv[1] = 1;
  rep(i, 2, n) { inv[i] = -(MOD / i) * inv[MOD % i] % MOD; }
  rep(i, 1, n) {
      b[i] = b[i - 1] * (k - 1 + i) % MOD * inv[i] % MOD;
      if(b[i]<0) b[i] += MOD;
  }
  conv(a, b, c, n);
  rep(i, 0, n) {
     printf("%lld\n", c[i]);
  }
}
int main() {
  int t = 1, ca = 0;

  while (t--) {

    solve();
  }
  return 0;
}
```
---

## 线性递推式

```c++
#include <bits/stdc++.h>
using namespace std;
#define rep(i,a,n) for (int i=a;i<n;i++)
#define per(i,a,n) for (int i=n-1;i>=a;i--)
#define pb push_back
#define mp make_pair
#define all(x) (x).begin(),(x).end()
#define fi first
#define se second
#define SZ(x) ((int)(x).size())
typedef vector<int> VI;
typedef long long ll;
typedef pair<int,int> PII;
const ll mod=1000000007;
ll powmod(ll a,ll b) {ll res=1;a%=mod; assert(b>=0); for(;b;b>>=1){if(b&1)res=res*a%mod;a=a*a%mod;}return res;}
// head

int _,n;
namespace linear_seq {
    const int N=10010;
    ll res[N],base[N],_c[N],_md[N];

    vector<int> Md;
    void mul(ll *a,ll *b,int k) {
        rep(i,0,k+k) _c[i]=0
        rep(i,0,k) if (a[i]) rep(j,0,k) _c[i+j]=(_c[i+j]+a[i]*b[j])%mod;
        for (int i=k+k-1;i>=k;i--) if (_c[i])
            rep(j,0,SZ(Md)) _c[i-k+Md[j]]=(_c[i-k+Md[j]]-_c[i]*_md[Md[j]])%mod;
        rep(i,0,k) a[i]=_c[i];
    }
    int solve(ll n,VI a,VI b) { // a 系数 b 初值 b[n+1]=a[0]*b[n]+...
//        printf("%d\n",SZ(b));
        ll ans=0,pnt=0;
        int k=SZ(a);
        assert(SZ(a)==SZ(b));
        rep(i,0,k) _md[k-1-i]=-a[i];_md[k]=1;
        Md.clear();
        rep(i,0,k) if (_md[i]!=0) Md.push_back(i);
        rep(i,0,k) res[i]=base[i]=0;
        res[0]=1;
        while ((1ll<<pnt)<=n) pnt++;
        for (int p=pnt;p>=0;p--) {
            mul(res,res,k);
            if ((n>>p)&1) {
                for (int i=k-1;i>=0;i--) res[i+1]=res[i];res[0]=0;
                rep(j,0,SZ(Md)) res[Md[j]]=(res[Md[j]]-res[k]*_md[Md[j]])%mod;
            }
        }
        rep(i,0,k) ans=(ans+res[i]*b[i])%mod;
        if (ans<0) ans+=mod;
        return ans;
    }
    VI BM(VI s) {
        VI C(1,1),B(1,1);
        int L=0,m=1,b=1;
        rep(n,0,SZ(s)) {
            ll d=0;
            rep(i,0,L+1) d=(d+(ll)C[i]*s[n-i])%mod;
            if (d==0) ++m;
            else if (2*L<=n) {
                VI T=C;
                ll c=mod-d*powmod(b,mod-2)%mod;
                while (SZ(C)<SZ(B)+m) C.pb(0);
                rep(i,0,SZ(B)) C[i+m]=(C[i+m]+c*B[i])%mod;
                L=n+1-L; B=T; b=d; m=1;
            } else {
                ll c=mod-d*powmod(b,mod-2)%mod;
                while (SZ(C)<SZ(B)+m) C.pb(0);
                rep(i,0,SZ(B)) C[i+m]=(C[i+m]+c*B[i])%mod;
                ++m;
            }
        }
        return C;
    }
    int gao(VI a,ll n) {
        VI c=BM(a);
        c.erase(c.begin());
        rep(i,0,SZ(c)) c[i]=(mod-c[i])%mod;
        return solve(n,c,VI(a.begin(),a.begin()+SZ(c)));
    }
};

int main() {
    for (scanf("%d",&_);_;_--) {
        scanf("%d",&n);
        printf("%d\n",linear_seq::gao(VI{2,24,96,416,1536,5504,18944,64000,212992,702464},n-1));
    }
}
```

---
## 高精度
```c++
#include<cstdio>
#include<cstring>
#include<iostream>
using std::max;
const int MAXL = 1000;
const int MAXB = 10000;
struct Bigint{
    int d[MAXL], l;
    void print() {
        printf("%d", d[l-1]);
        for(int i=l-2; i>-1; --i)
            printf("%04d", d[i]);
        printf("\n");
    }
    Bigint operator = (char *s) {
        memset(d, 0, sizeof(d));
        l = 0;
        for(int i=strlen(s)-1; i>-1; ++l)
            for(int j=1; j<MAXB&&i>-1; j *= 10, --i)
                d[l] += j*(s[i]-'0');
        return *this;
    }
    Bigint operator = (int b) {
        memset(d, 0, sizeof(d));
        l = 0;
        while(b) {
            d[l++] = b%MAXB;
            b /= MAXB;
        }
        if(!l) l = 1;
        return *this;
    }
    friend Bigint operator + (
            const Bigint& a,
            const Bigint& b) {
        Bigint c;
        c.l = max(a.l,b.l);
        int x = 0;
        for(int i=0; i<c.l; ++i) {
            x += a.d[i]+b.d[i];
            c.d[i] = x%MAXB;
            x /= MAXB;
        }
        if(x) c.d[c.l++] = x;
        return c;
    }
```

---

```c++
    friend Bigint operator * (
            const Bigint& a,
            const Bigint& b) {
        Bigint c;
        c.l = a.l+b.l;
        for(int i=0; i<a.l; ++i)
            for(int j=0; j<b.l; ++j)
                c.d[i+j] += a.d[i]*b.d[j];
        for(int i=0; i<c.l; ++i)
            if(c.d[i]>MAXB) {
                c.d[i+1] += c.d[i]/MAXB;
                c.d[i] %= MAXB;
            }
        while(c.l>1&&!c.d[c.l-1]) --c.l;
        return c;
    }
    friend Bigint operator * (
            const Bigint& a,
            const int& b) {
        Bigint c;
        c.l = a.l;
        int x = 0;
        for(int i=0; i<c.l; ++i) {
            x += b*a.d[i];
            c.d[i] = x%MAXB;
            x /= MAXB;
        }
        while(x) {
            c.d[c.l++] = x%MAXB;
            x /= MAXB;
        }
        return c;
    }
    friend Bigint operator - (
            const Bigint& a,
            const Bigint& b) {
        Bigint c;
        c.l = a.l;
        for(int i=0; i<c.l; ++i) {
            c.d[i] += a.d[i] - b.d[i];
            if(c.d[i]<0) {
                c.d[i] += MAXB;
                --c.d[i+1];
            }
        }
        while(c.l>1&&!c.d[c.l-1]) --c.l;
        return c;
    }
```

---

```c++
    friend Bigint operator / (
            const Bigint& a,
            const int& b) {
        Bigint c;
        c.l = a.l;
        int r = 0;
        for(int i=c.l-1; i>-1; --i) {
            r = r*MAXB+a.d[i];
            c.d[i] = r/b;
            r %= b;
        }
        while(c.l>1&&!c.d[c.l-1]) --c.l;
    }
    int cmp(const Bigint& b) {
        if(l<b.l) return -1;
        if(l>b.l) return 1;
        for(int i=l-1; i>-1; --i)
            if(d[i]-b.d[i]) {
                if(d[i]<b.d[i]) return -1;
                else return 1;
            }
        return 0;
    }
    friend Bigint operator / (
            const Bigint& a,
            const Bigint& b) {
        Bigint c, temp;
        c.l = a.l-b.l+1;
        for(int i=c.l-1; i>-1; --i) {
            if(i==c.l-1) {
                memcpy(temp.d, a.d+i, sizeof(int)*b.l);
                temp.l = b.l;
            }else{
                temp = temp*MAXB;
                temp.d[0] = a.d[i];
            }
            if(temp.cmp(b)==-1) c.d[i] = 0;
            else {
                int l=0, r=MAXB;
                while(l+1!=r) {
                    int mid = (l+r)/2;
                    int res = temp.cmp(b*mid);
                    if(res==0) {
                        l = mid;
                        break;
                    }
                    if(res==1) l = mid;
                    else r = mid;
                }
                c.d[i] = l;
                temp = temp-b*l;
            }
        }
        while(c.l>1&&!c.d[c.l-1]) --c.l;
        return c;
    }
}A, B;
```

---

# 图论
## LCA tarjian
```c++
#include<iostream>
#include<vector>
using namespace std;
const int MAXN = 4e4+1e3;
struct QNode{
    int u, v, lca;
}qe[1000];
struct Edge{
    int v, d;
};
int n, m, fa[MAXN], dis[MAXN];
bool vis[MAXN];
vector <Edge> G[MAXN];
vector <int> Query[MAXN];
void init(int n) {
    for(int i=1; i<=n; ++i) {
        fa[i] = i;
        G[i].clear();
        Query[i].clear();
        vis[i] = false;
        dis[i] = 0;
    }
}
```

---

```c++
int getfa(int x) {
    if(x==fa[x]) return x;
    else return fa[x]=getfa(fa[x]);
}
void dfs(int u, int pa) {
    for(int i=0; i<G[u].size(); ++i)
        if(G[u][i].v!=pa){
            int v = G[u][i].v;
            dis[v] = dis[u]+G[u][i].d;
            dfs(v, u);
            fa[v] = u;
        }
    vis[u] = true;
    for(int i=0; i<Query[u].size(); ++i) {
        int j = Query[u][i], v;
        if(qe[j].u==u) v = qe[j].v;
        else v = qe[j].u;
        if(vis[v]) qe[j].lca=getfa(v);
    }
}
int main() {
    ios::sync_with_stdio(false);
    int T;cin >> T;
    while(T--) {
        cin>> n>> m;
        init(n);
        for(int i=1; i<n; ++i) {
            int u, v, d;
            cin>> u>> v>> d;
            G[u].push_back((Edge){v,d});
            G[v].push_back((Edge){u,d});
        }
        for(int i=0; i<m; ++i) {
            cin>> qe[i].u>> qe[i].v;
            Query[qe[i].u].push_back(i);
            Query[qe[i].v].push_back(i);
        }
        dfs(1,0);
        for(int i=0; i<m; ++i)
            cout<< (dis[qe[i].u]-dis[qe[i].lca])+
                    (dis[qe[i].v]-dis[qe[i].lca])<<'\n';
    }
    return 0;
}
```

---

# <big>字符串</big>
## 后缀数组
```c++
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
const int MAXN = 1e5+1e3;
typedef long long ll;
int sa[MAXN];
int r[MAXN], height[MAXN];
bool cmp(int *r, int a, int b, int l) {
    return r[a]==r[b]&&r[a+l]==r[b+l];
}
void DA(char* s, int n) {
    static int wa[MAXN], wb[MAXN], cnt[MAXN];
    int *x=wa, *y=wb;
    int m=128, i, j, p;
    for(i=0; i<m; ++i) cnt[i]=0;
    for(i=0; i<n; ++i) ++cnt[x[i]=s[i]];
    for(i=1; i<m; ++i) cnt[i]+=cnt[i-1];
    for(i=n-1; i>-1; --i) sa[--cnt[x[i]]] = i;
    for(j=1; j<n; m=p,j*=2) {
        for(p=0,i=n-j; i<n; ++i)
            y[p++] = i;
        for(i=0; i<n; ++i) if(sa[i]>=j)
            y[p++] = sa[i]-j;
        for(i=0; i<m; ++i) cnt[i] = 0;
        for(i=0; i<n; ++i) ++cnt[x[i]];
        for(i=1; i<m; ++i) cnt[i]+=cnt[i-1];
        for(i=n-1; i>-1; --i)
            sa[--cnt[x[y[i]]]] = y[i];
        swap(x,y);
        p = 1;x[sa[0]]=0;
        for(int i=1; i<n; ++i)
            x[sa[i]] = cmp(y,sa[i-1], sa[i], j)?p-1:p++;
        if(p==n) break;
    }
}
```

---

```c++
void calheight(char *s, int n) {
    int i, j, k=0;
    for(i=0; i<n; ++i)
        r[sa[i]] = i;
    for(i=0; i<n; ++i) {
        if(k) k--;
        j = sa[r[i]-1];
        while(s[i+k]==s[j+k]) k++;
        height[r[i]] = k;
    }
}
```
## 最长回文子串
```c++
#include<cstdio>
#include<cstring>
const int MAXN = 1e5+2e4;
void manacher(char *s) {
    int id=0, mx=0;
    for(int i=1; s[i]; ++i) {
        if(mx>i) p[i] = min(p[id*2-i],mx-i);
        else p[i] = 1;
        while(s[i+p[i]]==s[i-p[i]]) p[i]++;
        if(i+p[i]>mx) {
            id = i;
            mx = i+p[i];
        }
    }
}
void solve(char *s) {
    static char t[MAXN*2];
    int sz = 0;
    t[sz++] = '$';
    t[sz++] = '#';
    for(int i=0; s[i]; ++i) {
        t[sz++] = s[i];
        t[sz++] = '#';
    }
    t[sz] = '\0';
    manacher(t);
}
```
## ShiftOr
```c++
#include<cstdio>
#include<cstring>
#include<bitset>
using namespace std;
const int MAXN = 5e6+1000;
const int MAXM = 128;
bitset<MAXM> B[10], D;
char s[MAXN], t[MAXM];
int main() {
	int i, ret;
	scanf("%s", s);
	scanf("%s", t);
	int len = strlen(t);
	D.set();
	for(i=0; i<10; ++i)
		B[i].set();
	for(i=0; i<len; ++i)
		B[t[i]-'0'][i] = 0;
	for(i=0; s+i<p; ++i) {
		D = D<<1|B[s[i]-'0'];
		if(!D[len-1]) {
			printf("%d\n", i-len+2);
			return 0;
		}
	}
	return 0;
}
```
## LCP
```c++
#include<cstdio>
#include<cstring>
typedef unsigned long long uint;
const int MAXN = 1e5+1e3;
int n;
char s[MAXN];
uint h[MAXN], x[MAXN];
void initia() {
	for(int i=n-1; i>-1; --i)
		h[i] = h[i+1]*3+s[i];
	x[0] = 1;
	for(int i=1; i<=n; ++i)
		x[i] = x[i-1]*3;
}
bool check(int p, int q, int l) {
	return h[p]-h[p+l]*x[l]==h[q]-h[q+l]*x[l];
}
int max(int a, int b) {
	return a>b?a:b;
}
int lcp(int p, int l) {
	int q = p+l;
	int ll = 0, rr = n-q;
	while(ll!=rr) {
		int m = (ll+rr)/2+1;
		if(check(p,q,m)) ll=m;
		else rr=m-1;
	}
	return ll;
}
int main() {
	scanf("%s", s);
	n = strlen(s);
	initia();
	int res = 0;
	for(int l=1; l<=n; ++l)
		for(int i=0; i+l<=n; i+=l) {
			int len = lcp(i, l);
			res = max(res, len/l+1);
			res = max(res,lcp(i-l+len%l, l)/l+1);
		}
	printf("%d\n", res);
	return 0;
}
```


## SAM
```c++
const int MAXN = 2e5+1e3;
struct state {
    int len, pre, ch[26];
    ll cnt;
};
struct SAM {
    int sz, last;
    state st[MAXN];
    state& operator[] (int x) {
        return st[x];
    }
    SAM() {
        sz=1, last=0;
        st[0].len=0, st[0].pre=-1;
    }
    void add(int c) {
        int cur=sz++, p;
        st[cur].len=st[last].len+1;
        for(p=last; p!=-1&&!st[p].ch[c]; p=st[p].pre)
            st[p].ch[c]=cur;
        if(p==-1) st[cur].pre=0;
        else {
            int q=st[p].ch[c];
            if(st[q].len==st[p].len+1)
                st[cur].pre=q;
            else {
                int clone=sz++;
                st[clone]=st[q];
                st[clone].len=st[p].len+1;
                st[cur].pre=st[q].pre=clone;
                for(; p!=-1&&st[p].ch[c]==q; p=st[p].pre)
                    st[p].ch[c]=clone;
            }
        }
        last=cur;
    }
    int find(char *t) {//查询lcs
        int now=0, l=0, ans=0;
        for(char* p=t; *p; ++p) {
            while(now&&!st[now].ch[*p-'a']) {
                now=st[now].pre;
                l=st[now].len;
            }
            if(st[now].ch[*p-'a']) {
                ++l;
                now=st[now].ch[*p-'a'];
            }
            ans=max(l, ans);
        }
        return ans;
    }
} sam;
char s[MAXN];
ll dfs(int u) {
    ll& ret=sam[u].cnt;
    if(ret) return ret;
    ret = 1;
    rep(i, 0, 26)
    if(sam[u].ch[i])
        ret += dfs(sam[u].ch[i]);
    return ret;
}
void printstr(int now, int k) {//输出第K大
    while(--k) {
        rep(i, 0, 26) {
            int v=sam[now].ch[i];
            if(v) {
                if(k<=sam[v].cnt) {
                    putchar(i+'a');
                    now=v;
                    break;
                }
                else k-= sam[v].cnt;
            }
        }
    }
    puts("");
}
int main() {
    sf("%s", s);
    for(char* p=s; *p; ++p)
        sam.add(*p-'a');
    dfs(0);
    int q;
    sf("%d", &q);
    while(q--) {
        int k;
        sf("%d", &k);
        printstr(0, k+1);
    }
    return 0;
}
```

## 最小表示法
```c++
#include<cstdio>
#include<cstring>
const int MAXN = 2048;
char s[MAXN];
int l;
int min(int a, int b) {
    return a<b?a:b;
}
void print(int st) {
    for(int i=0; i<l; ++i)
        putchar(s[i+st]);
    puts("");
}
int main() {
    scanf("%s", s);
    l = strlen(s);
    for(int i=0; i<l; ++i)
        s[i+l] = s[i];
    s[l<<1] = 0;
    if(l==1) {
        print(0);
        return 0;
    }
    int i, j;
    for(i=0, j=1; i<l&&j<l; ) {
        int k=0;
        while(k<l&&s[i+k]==s[j+k]) ++k;
        if(k==l) break;
        if(s[i+k]>s[j+k]) i=i+k+1;
        else j=j+k+1;
        if(i==j) j+=1;
    }
    print(min(i,j));
    return 0;
}
```

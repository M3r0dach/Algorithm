<h1 align = "center">目录</h1>

1. DP
    - 斜率优化的DP
    - 大背包
1. 计算几何  
    - 凸包
    - 动态凸包
    - 旋转卡壳法
    - 矩形面积求并
1. 数论
    - 离散对数
    - 欧拉函数
    - 素数测试
    - 素数个数
    - Lucas定理
1. 区间问题
    - 莫队算法
    - ST表
    - fzu2224
    - Treap
    - 树上莫队cot2
1. 数学
	- FFT
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

#DP
##斜率优化的DP

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

##大背包
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

#计算几何
##凸包
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

##动态凸包
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

##旋转卡壳法
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
##矩形面积求并
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

* * *

#数论
##离散对数
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
##欧拉函数
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
##素数测试
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
##素数个数
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

##Lucas定理
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

#区间问题
##莫队算法
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
##ST表
```c++
#include<cstdio>
#include<iostream>
using namespace std;
const int MAXN = 200000;
const int MAXL = 20;
int st[MAXN][MAXL];
int n, m, l, r, len;
int main() {
    int i, j, k;
    scanf("%d%d", &n, &m);
    for(i = 1; i < n+1; ++i)
        scanf("%d", &st[i][0]);
    for(j = 1; 1<<j <= n; ++j)
        for(i = 1; i+(1<<j) < n+2; ++i)
            st[i][j] = min(st[i][j-1],st[i+(1<<(j-1))][j-1]);
    for(i = 0; i < m; ++i) {
        scanf("%d%d", &l, &r);
        len = r-l+1, k = 0;
        while((1<<(k+1)) < len) ++k;
        printf("%d ", min(st[l][k],st[r+1-(1<<k)][k]));
    }
    return 0;
}
```
##fzu2224
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
##Treap
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
##树上莫队cot2
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
#数学
##FFT
```c++
#include<bits/stdc++.h>
using namespace std;
const int MaxN = 262144;
const double PI = 3.14159265358979;
complex<double> fft_pow[MaxN];
void fft(
       const complex<double> *a,
       const int &n, const int &step,
       complex<double> *out,
       const int type)
{
   if (n == 1) {
       out[0] = a[0];
       return;
   }
   int m = n >> 1;
   fft(a, m, step << 1, out, type);
   fft(a + step, m, step << 1, out + m, type);
   for (int i = 0; i < m; i++) {
       complex<double> even = out[i];
       complex<double> odd = out[i + m];
       if (type == 1) odd *= fft_pow[i * step];
       else odd /= fft_pow[i * step];
       out[i] = even + odd;
       out[i + m] = even - odd;
   }
}
void fft(
       const complex<double> *a, const complex<double> *b,
       int n,
       complex<double> *res)
{
   int tn;
   for (tn = 1; tn < n + n; tn <<= 1);
   n = tn;
   for (int i = 0; i < n; i++)
       fft_pow[i] = exp(complex<double>(0.0, -2 * PI * i / n));
   static complex<double> ta[MaxN];
   static complex<double> tb[MaxN];
   fft(a, n, 1, ta, 1);
   fft(b, n, 1, tb, 1);
   for (int i = 0; i < n; i++)
       ta[i] *= tb[i];
   fft(ta, n, 1, res, -1);
   for (int i = 0; i < n; i++)
       res[i] /= n;
}
```

---

##高精度
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
                c.d[i+j] += a.d[i]*b.d[i];
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

#图论
##LCA tarjian
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

#<big>字符串</big>
##后缀数组
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
##最长回文子串
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
##ShiftOr
```
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
##LCP
```
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

##SAM
```
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
```
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

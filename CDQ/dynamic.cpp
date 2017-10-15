#include<cstdio>
#include<cstring>
#include<algorithm>
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = 1e5+1e3;
struct Act{
    int t, pos, v;
}act[MAXN];
int n, m, bit[MAXN], sz;
int a[MAXN], id[MAXN];
ll f[MAXN], ans[MAXN];
bool used[MAXN];
void add(int x, int v) {
    while(x<n+1) {
        bit[x] += v;
        x += x&(-x);
    }
}
int sum(int x) {
    int ret=0;
    while(x) {
        ret += bit[x];
        x -= x&(-x);
    }
    return ret;
}

bool cmp(const Act& a,const Act& b) {
    return a.pos<b.pos;
}
void cdq(int l, int r) {
    if(l==r) return;
    int mid=(l+r)/2;
    cdq(l,mid);cdq(mid+1,r);
    sort(act+l, act+r+1, cmp);
    sz=0;
    rep(i,l,r+1) {
        if(act[i].t<=mid) {
            add(act[i].v, 1);
            ++sz;
        }
        else f[act[i].t] += sz-sum(act[i].v);
    }
    rep(i,l,r+1)
        if(act[i].t<=mid)
            add(act[i].v, -1);

    drep(i, r, l-1) {
        if(act[i].t<=mid) {
            add(act[i].v, 1);
            ++sz;
        }
        else f[act[i].t] += sum(act[i].v);
    }
    rep(i,l,r+1)
        if(act[i].t<=mid)
            add(act[i].v, -1);
}
int main() {
    sf("%d%d", &n, &m);
    rep(i,1,n+1) {
        sf("%d", a+i);
        id[a[i]] = i;
    }
    int k=n;
    rep(i,1,m+1) {
        int v; sf("%d", &v);
        act[k]=(Act){k, id[v], v};
        used[id[v]] = true;
        --k;
    }
    rep(i,1,n+1) if(!used[i]) {
        act[k]=(Act){k, i, a[i]};
        --k;
    }
    cdq(1,n);
    rep(i,1,n+1)
        ans[i] = ans[i-1]+f[i];
    drep(i,n,n-m)
        pf("%lld\n", ans[i]);
    return 0;
}

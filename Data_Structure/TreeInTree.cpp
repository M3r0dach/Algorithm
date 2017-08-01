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

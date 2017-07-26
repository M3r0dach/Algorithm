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

//HDU6040
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
const int MAXN = 1e7+1e3;
const int MAXM = 128;
typedef pair<int,int> pii;
unsigned n, m, A, B, C;
unsigned x, y, z;
unsigned a[MAXN], ans[MAXM];
pii query[MAXM];
unsigned rng61() {
    unsigned t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;
    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;
    return z;
}
int saperate(unsigned* a, unsigned sz) {
    unsigned pivot=a[0], i=0,l=0, r=sz-1;
    while(l<r) {
        while(l<r&&a[r]>a[i]) --r;
        swap(a[r],a[i]);
        i=r;
        if(l==r) break;
        while(l<r&&a[l]<=a[i]) ++l;
        swap(a[l],a[i]);
        i=l;
        if(l==r) break;
    }
    return i;
}
void work(unsigned* a, unsigned asz, pii* q, unsigned qsz) {
    if(asz<=0||qsz<=0) return;
    if(asz<=3) {
        for(unsigned i=0; i<asz; ++i)
            for(unsigned j=i+1; j<asz; ++j)
                if(a[i]>a[j]) swap(a[i],a[j]);
        for(unsigned i=0; i<qsz; ++i)
            ans[q[i].second] = a[q[i].first];
        return;
    }
    int lsz = saperate(a,asz);
    int m1=0, m2=qsz;
    for(int i=0; i<qsz; ++i) {
        if(q[i].first<lsz) {
            m1=i+1;
        } 
        else if(q[i].first==lsz) {
            ans[q[i].second] = a[lsz];
        }
        else {
            m2=min(m2,i);
            q[i].first-=lsz+1;
        }
    }
    work(a,lsz, q, m1);
    work(a+lsz+1,asz-lsz-1,q+m2, qsz-m2);
}
int main() {
    int ca=0;
    while(~scanf("%u%u%u%u%u", &n, &m, &A, &B, &C)) {
        x = A, y = B, z = C;
        for(int i=0; i<m; ++i) {
            int b; scanf("%d", &b);
            query[i] = make_pair(b,i);
        }
        for(unsigned i=0; i<n; ++i)
            a[i] = rng61();
        sort(query, query+m);
        work(a,n,query,m);
        printf("Case #%d:", ++ca);
        for(unsigned i=0; i<m; ++i)
            printf(" %u", ans[i]);
        puts("");
    }
    return 0;
}


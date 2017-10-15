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

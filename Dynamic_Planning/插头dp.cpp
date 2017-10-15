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

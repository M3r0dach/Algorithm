//https://vjudge.net/contest/169617#problem/J
#include<cstdio>
#include<cstring>
#include<cmath>
#include<algorithm>
using namespace std;
const int MAXN = 56;
typedef long long ll;
bool g[MAXN][MAXN];
ll a[MAXN][MAXN];
ll eliminate(int n) {
    ll res=1;
    for(int i=0; i<n; ++i) {
        for(int j=i+1; j<n; ++j) {
            int x=i, y=j;
            while(a[y][i]) {
                ll t = a[x][i]/a[y][i];
                for(int k=0; k<n; ++k)
                    a[x][k] -= a[y][k]*t;
                swap(x,y);
            }
            if(x!=i) {
                for(int k=0; k<n; ++k)
                    swap(a[x][k], a[i][k]);
            }
        }
        if(a[i][i]==0) return 0;
        res *= a[i][i];
    } 
    if(res<0) res*=-1;
    return res;
}
int main() {
    int n, m, t;
    scanf("%d", &t);
    while(t--) {
        scanf("%d%d", &n, &m);
        memset(g, false, sizeof(g));
        memset(a, 0, sizeof(a));
        while(m--) {
            int u, v;
            scanf("%d%d", &u, &v);
            --u, --v;
            g[u][v] = g[v][u] = true;
        }
        for(int i=0; i<n; ++i)
            for(int j=i+1; j<n; ++j)
                if(g[i][j]) {
                    ++a[i][i], ++a[j][j];
                    a[i][j] = a[j][i] = -1;
                }
        printf("%lld\n", eliminate(n-1));
    }
    return 0;
}

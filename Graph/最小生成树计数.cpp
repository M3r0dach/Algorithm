//https://vjudge.net/problem/HYSBZ-1016
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
using namespace std;
const int MAXN = 128;
const int MAXM = 1024;
const int MOD = 31011;

struct Edge{
    int u, v, d;
    bool operator < (const Edge& b) const
    {
        return d<b.d;
    }
}e[MAXM];
int unit[MAXN], fa[MAXN];
int n, m;
bool vis[MAXN];
int A[MAXN][MAXN], B[MAXN][MAXN];
vector<pair<int,int> > g[MAXN];
vector<int> v;

int findroot(int x, int *p) {
    if(x==p[x]) return x;
    return p[x]=findroot(p[x], p);
}
int det(int A[][MAXN], int n) {
    int ret=1, sgn=0;
    for(int i=0; i<n; ++i) {
        for (int j = i+1; j <n; ++j) {
            int x=i, y=j;
            while (A[y][i]) {
                int t=A[x][i]/A[y][i]%MOD;
                for(int k=i; k<n; ++k)
                    A[x][k] = (A[x][k]-A[y][k]%MOD*t)%MOD;
                swap(x,y);
            }
            if(x!=i) {
                for(int k=i; k<=n; ++k)
                    swap(A[x][k], A[y][k]);
                ++sgn;
            }
        }
        if(A[i][i]==0) return 0;
        else ret=ret*A[i][i]%MOD;
    }
    if(sgn&1) ret=-ret;
    return ret;
}
int make_graph(int c) {
    memset(vis, false, sizeof(vis));
    memset(A, 0, sizeof(A));
    memset(B, 0, sizeof(B));
    for(int i=0; i<g[c].size(); ++i) {
        int u=g[c][i].first;
        int v=g[c][i].second;
        vis[u]=vis[v]=true;
        ++B[u][v], ++B[v][u];
    }
    v.clear();
    for(int i=1; i<=n; ++i)
        if(vis[i]) v.push_back(i);
    for(int i=0; i<v.size(); ++i)
        for(int j=i+1; j<v.size(); ++j)
        {
            int a=v[i], b=v[j];
            if(B[a][b]) {
                A[i][j] -= B[a][b];
                A[j][i] -= B[a][b];
                A[i][i] += B[a][b];
                A[j][j] += B[a][b];
            }
        }
    return det(A, v.size()-1);
}
int main() {
    scanf("%d%d", &n, &m);
    for(int i=0; i<m; ++i)
        scanf("%d%d%d", &e[i].u, &e[i].v, &e[i].d);
    sort(e, e+m);
    for(int i=1; i<=n; ++i)
        fa[i] = unit[i] = i;
    int res=1;
    for(int i=0; i<m; ) {
        int j=i+1;
        while (j<m&&e[j].d==e[i].d) ++j;
        for(int k=i; k<j; ++k) {
            int u=findroot(e[k].u, fa);
            int v=findroot(e[k].v, fa);
            if(u!=v) fa[u] = v;
        }
        for(int k=1; k<=n; ++k)
            g[k].clear();
        for(int k=i; k<j; ++k) {
            int u=findroot(e[k].u, unit);
            int v=findroot(e[k].v, unit);
            int c=findroot(e[k].u, fa);
            if(u!=v) g[c].push_back(make_pair(u,v));
        }
        for(int k=1; k<=n; ++k)
            if(g[k].size())
                res = res*make_graph(k)%MOD;
        for(int k=1; k<=n; ++k)
            unit[k] =fa[k];
        i=j;
        if(res==0) break;
    }
    int cnt=0;
    for(int i=1; i<=n; ++i)
        if(fa[i]==i) ++cnt;
    if(cnt>1) res = 0;
    if(res<0) res += MOD;
    printf("%d\n", res);
    return 0;
}

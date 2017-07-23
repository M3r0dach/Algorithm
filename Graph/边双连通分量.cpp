//POJ3352
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<vector>
#include<cmath>
#define fi first
#define se second
using namespace std;
const int MAXN = 1024;
struct Edge{
    int u, v;
    int pal(int x) {
        return u+v-x;
    }
}e[MAXN];
int dfs_clk;
int pre[MAXN], low[MAXN], rk[MAXN];
int n, r;
int fa[MAXN], deg[MAXN];
vector<int> G[MAXN];

int getfa(int x) {
    if(x==fa[x]) return x;
    return fa[x]=getfa(fa[x]);
}

void dfs(int u, int pe) {
    pre[u]=low[u]=++dfs_clk;
    rk[dfs_clk] = u;
    for(int i=0; i<G[u].size(); ++i) {
        if(G[u][i]==pe) continue;
        int v=e[G[u][i]].pal(u);
        if(!pre[v]) {
            dfs(v, G[u][i]);
            low[u] = min(low[u], low[v]);
        }
        else low[u] = min(low[u], pre[v]);
    }
}
int main() {
    scanf("%d%d", &n, &r);
    dfs_clk = 0;
    memset(pre, 0, sizeof(pre));
    memset(deg, 0, sizeof(deg));
    for(int i=1; i<=n; ++i) {
        G[i].clear();
        fa[i] = i;
    }
    for(int i=0; i<r; ++i) {
        scanf("%d%d", &e[i].u, &e[i].v);
        G[e[i].u].push_back(i);
        G[e[i].v].push_back(i);
    }
    dfs(1, -1);
    for(int i=1; i<=n; ++i) {
        int u=getfa(rk[low[i]]), v=getfa(i);
        if(u!=v) fa[u] = v;
    }
    for(int i=0; i<r; ++i) {
        int u=getfa(e[i].u);
        int v=getfa(e[i].v);
        if(u==v) continue;
        ++deg[u], ++deg[v];
    }
    int cnt=0;
    for(int i=1; i<=n; ++i)
        if(deg[i]==1) ++cnt;
    printf("%d\n", (cnt+1)/2);
    return 0;
}

//POJ 1523
#include<cstdio>
#include<cstring>
#include<vector>
#include<algorithm>
using namespace std;
const int MAXN = 1024;
struct Edge {
    int u, v;
    int pal(int x) {
        return u+v-x;
    }
};
int pre[MAXN], low[MAXN], dfs_clk;
int cnt[MAXN];
vector<int> G[MAXN], spf;
vector<Edge> edges;
void init() {
    edges.clear();
    spf.clear();
    for(int i=0; i<MAXN; ++i)
        G[i].clear();
    dfs_clk=0;
    memset(pre, 0, sizeof(pre));
    memset(cnt, 0, sizeof(cnt));
}
void dfs(int u, int pe) {
    pre[u]=low[u]=++dfs_clk;
    int chsz=0;
    int cut=0;
    for(int i=0; i<G[u].size(); ++i) {
        if(G[u][i]==pe) continue;
        int v=edges[G[u][i]].pal(u);
        if(!pre[v]) {
            ++chsz;
            dfs(v, G[u][i]);
            if(low[v]>=pre[u]) ++cut;
            low[u] = min(low[u], low[v]);
        }
        else low[u] = min(low[u], pre[v]);
    }
    if(pe==-1&&chsz<2) cut=0;
    if(cut) {
        spf.push_back(u);
        cnt[u] = cut;
        if(u!=1) ++cnt[u];
    }
}
void solve() {
    static int cas=0;
    for(int i=1; i<MAXN; ++i)
        if(!pre[i]) dfs(i, -1);
    printf("Network #%d\n", ++cas);
    if(spf.size()) {
        sort(spf.begin(), spf.end());
        for(int i=0; i<spf.size(); ++i)
            printf("  SPF node %d leaves %d subnets\n", spf[i], cnt[spf[i]]);
    }
    else puts("  No SPF nodes");
    puts("");
}
int main() {
    int u, v;
    while(~scanf("%d", &u)) {
        if(u==0) {
            if(edges.size()==0)
                break;
            solve();
            init();
            continue;
        }
        scanf("%d", &v);
        int m = edges.size();
        edges.push_back((Edge) {u,v});
        G[u].push_back(m);
        G[v].push_back(m);
    }
    return 0;
}

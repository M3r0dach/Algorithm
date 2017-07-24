#include<cstdio>
#include<cstring>
#include<stack>
#include<vector>
using namespace std;
const int MAXN = 1e4+1e3;
int dfs_clk, pre[MAXN], low[MAXN];
int scco[MAXN], scco_cnt;
vector<int> G[MAXN];
stack<int> S;
int n, m;
void dfs(int u) {
    pre[u]=low[u]=++dfs_clk;
    S.push(u);
    for(int i=0; i<G[u].size(); ++i) {
        int v=G[u][i];
        if(!pre[v]) {
            dfs(v);
            low[u] = min(low[u], low[v]);
        }
        else if(!scco[v])
            low[u] = min(low[u], pre[v]);
    }
    if(low[u]==pre[u]) {
        ++scco_cnt;
        while(S.top()!=u) {
            int v=S.top();
            S.pop();
            scco[v] = scco_cnt;
        }
        scco[u] = scco_cnt;
        S.pop();
    }
}
int main() {
    while(~scanf("%d%d", &n, &m)) {
        if(!n&&!m) break;
        for(int i=0; i<=n; ++i)
            G[i].clear();
        while(!S.empty()) S.pop();
        memset(scco, 0, sizeof(scco));
        memset(pre, 0, sizeof(pre));
        dfs_clk = scco_cnt = 0;
        for(int i=0; i<m; ++i) {
            int u, v;
            scanf("%d%d", &u, &v);
            G[u].push_back(v);
        }
        for(int i=1; i<=n; ++i)
            if(!scco[i]) dfs(i);
        if(scco_cnt>1) puts("No");
        else puts("Yes");
    }
    return 0;
}

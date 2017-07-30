//HDU 6041
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<vector>
#include<stack>
#include<queue>
using namespace std;
const int MAXN = 1024;
struct Edge{
    int u, v, d;
    Edge(int uu=0, int vv=0, int dd=0)
        :u(uu), v(vv), d(dd){}
    int pal(int x) {
        return u^v^x;
    }
}edges[MAXN*2];
vector<int> G[MAXN];
vector<int> res, W;
int low[MAXN], pre[MAXN];
int n, m, dfs_clk, cnt, K;
int sum;
stack<Edge> S;
void init() {
    for(int i=1; i<=n; ++i)
        G[i].clear();
    cnt=dfs_clk=0;
    memset(pre, 0, sizeof(pre));
    res.clear();
}
struct Opt{
    int w, x, y;
    Opt(int ww=0, int xx=0, int yy=0)
        :w(ww), x(xx), y(yy){}
    bool operator <(const Opt& b) const {
        if(w!=b.w) return w<b.w;
        if(x!=b.x) return x<b.x;
        if(y!=b.y) return y<b.y;
        return false;
    }
};
void myMerge(vector<int>& V, vector<int>& b) {
    priority_queue<Opt> pq;
    for(int i=0; i<b.size(); ++i)
        pq.push(Opt(V[0]+b[i], i, 0));
    W.clear();
    while(W.size()<K&&!pq.empty()) {
        auto e=pq.top();
        pq.pop();
        W.push_back(e.w);
        if(e.y+1<V.size()) 
            pq.push(Opt(V[e.y+1]+b[e.x],e.x,e.y+1));
    }
    V = W;
}
void dfs(int u, int pe) {
    pre[u]=low[u]=++dfs_clk;
    for(int i=0; i<G[u].size(); ++i) {
        if(G[u][i]==pe) continue;
        Edge& e=edges[G[u][i]];
        int v=e.pal(u);
        if(!pre[v]) {
            S.push(e);
            dfs(v, G[u][i]);
            low[u] = min(low[u], low[v]);
            if(low[v]>=pre[u]) {
                vector<int> vec;
                while(!S.empty()) {
                    Edge sel = S.top();
                    S.pop();
                    vec.push_back(sel.d);
                    if(sel.u==u&&sel.v==v||sel.u==v&&sel.v==u)
                        break;
                }
                if(vec.size()>1) myMerge(res, vec);
            }
        }
        else if(pre[v]<pre[u]){//important
            S.push(e);
            low[u] = min(low[u], pre[v]);
        }
    }
}
int main() {
    int ca=0;
    res.reserve(100100);
    while(~scanf("%d%d", &n, &m)) {
        init();
        sum = 0;
        for(int i=0; i<m; ++i) {
            int u, v, d;
            scanf("%d%d%d", &u, &v, &d);
            edges[i] = Edge(u,v,d);
            G[u].push_back(i);
            G[v].push_back(i);
            sum += d;
        }
        scanf("%d", &K);
        res.push_back(0);
        dfs(1,-1);
        int ans=0;
        for(int i=0; i<res.size(); ++i)
            ans += (i+1)*(sum-res[i]);
        printf("Case #%d: %u\n", ++ca, ans);
    }
    return 0;
}

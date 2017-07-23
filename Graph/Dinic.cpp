//HDU - 3572 
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<queue>
#include<vector>
using namespace std;
const int MAXN = 2048;
const int MAXD = 512;
const int INF = 0x3f3f3f3f;

struct Edge{
    int fr, to, cap;
};
vector<int> G[MAXN];
vector<Edge> edges;
void init() {
    edges.clear();
    for(int i=0; i<MAXN; ++i)
        G[i].clear();
}
void addedge(int fr, int to, int cap) {
    edges.push_back((Edge){fr, to, cap});
    edges.push_back((Edge){to, fr, 0});
    int m=edges.size();
    G[fr].push_back(m-2);
    G[to].push_back(m-1);
}
int d[MAXN], cur[MAXN];
int bfs(int s, int t) {
    memset(d, 0, sizeof(d));
    queue<int> Q;
    d[s]=1;
    Q.push(s);
    while(!Q.empty()) {
        int u=Q.front();
        Q.pop();
        for(int i=0; i<G[u].size(); ++i) {
            Edge& e=edges[G[u][i]];
            if(!d[e.to]&&e.cap>0) {
                d[e.to]=d[u]+1;
                Q.push(e.to);
            }
        }
    }
    return d[t];
}
int dfs(int u, int t, int a) {
    if(!a||u==t) return a;
    int flow=0;
    for(int& i=cur[u]; i<G[u].size(); ++i) {
        Edge& e=edges[G[u][i]];
        if(d[u]+1==d[e.to]) {
            int f=dfs(e.to, t, min(e.cap, a));
            if(f>0) {
                flow += f;
                e.cap -= f;
                edges[G[u][i]^1].cap += f;
                a -= f;
                if(!a) break;
            }
        }
    }
    return flow;
}
int Dinic(int s, int t) {
    int ret=0;
    while(bfs(s,t)) {
        memset(cur, 0, sizeof(cur));
        ret += dfs(s,t,INF);
    }
    return ret;
}
int main() {
    int t; scanf("%d", &t);
    int cas=0;
    while(t--) {
        int n, m;
        scanf("%d%d", &n, &m);
        int sum=0;
        init();
        for(int i=1; i<=n; ++i) {
            int p, s, e;
            scanf("%d%d%d", &p, &s, &e);
            sum += p;
            addedge(0,i,p);
            for(int j=s; j<=e; ++j)
                addedge(i, n+j, 1);
        }
        for(int i=1; i<=MAXD; ++i)
            addedge(n+i, n+MAXD+1, m);
        if(sum==Dinic(0, n+MAXD+1)) printf("Case %d: Yes\n\n", ++cas);
        else printf("Case %d: No\n\n", ++cas);
    }
    return 0;
}

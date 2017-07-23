#include<cstdio>
#include<vector>
#include<cstring>
#include<algorithm>
using namespace std;
const int MAXN = 1e4+1e3;
struct Edge{
    int to, l;
    Edge* next;
    Edge() {}
    Edge(int to, int l, Edge* next)
        :to(to), l(l), next(next){}
}e[MAXN*2], *head[MAXN];
vector<int> dep;
int n, k, size;
int sz[MAXN], f[MAXN], root;
bool done[MAXN];

void init() {
    memset(done, false, sizeof(done));
    memset(head, 0, sizeof(head));
    f[0] = n+1;
}
void getroot(int u, int fa) {
    sz[u]=1, f[u]=0;
    for(Edge* p=head[u]; p; p=p->next) {
        int v=p->to;
        if(done[v]||v==fa) continue;
        getroot(v, u);
        sz[u] += sz[v];
        f[u] = max(f[u], sz[v]);
    }
    f[u] = max(f[u], size-sz[u]);
    if(f[u]<f[root]) root=u;
}
void getdep(int u, int fa, int d) {
    dep.push_back(d);
    sz[u]=1;
    for(Edge* p=head[u]; p; p=p->next) {
        int v=p->to;
        if(done[v]||v==fa) continue;
        getdep(v, u, d+p->l);
        sz[u]+=sz[v];
    }
}
int cal(int node, int rh) {
    if(k<2*rh) return 0;
    dep.clear();
    getdep(node, -1, rh);
    sort(dep.begin(), dep.end());
    int r=dep.size()-1, l=0, ret=0;
    while(l<r) {
        if(dep[l]+dep[r]>k) --r;
        else ret += r-l++;
    }
    return ret;
}
int work(int node, int n) {
    int ret=0;
    size = n;
    root=0;
    getroot(node, -1);
    node=root;

    ret = cal(node, 0);
    done[node] = true;
    for(Edge* p=head[node]; p; p=p->next) {
        if(done[p->to]) continue;
        ret -= cal(p->to, p->l);
        ret += work(p->to, sz[p->to]);
    }
    return ret;
}
int main() {
    while(~scanf("%d%d", &n, &k)) {
        if(!n&&!k) break;
        init();
        for(int i=0; i<n-1; ++i) {
            int u, v, l;
            scanf("%d%d%d", &u, &v, &l);
            e[i<<1] = Edge(v,l,head[u]);
            head[u] = e+i*2;
            e[i<<1|1] = Edge(u,l,head[v]);
            head[v] = e+i*2+1;
        }
        printf("%d\n", work(1,n));
    }
    return 0;
}

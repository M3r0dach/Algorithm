#include<cstdio>
#include<cstring>
#include<algorithm>
#include<vector>
#include<cmath>
using namespace std;
const int MAXN = 7e4;
typedef pair<int,int> pii;
struct Edge{
    int to, next;
    Edge(int to=0, int next=0)
        :to(to), next(next){};
}e[MAXN];
int head[MAXN], esz;
int fp[MAXN];
struct Node{
    int top, fa, dep;
    int num, pos, son;
    int v;
}node[MAXN];
int n, pos;
pii tr[MAXN*15];

pii operator * (const pii& a, const pii& b) {
    return make_pair(max(a.first, b.first), a.second+b.second);
}

void init() {
    for(int i=1; i<=n; ++i)
        head[i]=node[i].son=-1;
    pos=1;
}
void dfs1(int u, int fa, int d) {
    node[u].fa=fa;
    node[u].dep=d;
    node[u].num=1;
    for(int p=head[u]; p!=-1; p=e[p].next) {
        int v=e[p].to;
        if(v==fa) continue;
        dfs1(v, u, d+1);
        node[u].num+=node[v].num;
        if(node[u].son==-1||
                node[node[u].son].num<node[v].num)
            node[u].son=v;
    }
}
void getpos(int u, int sp) {
    node[u].top=sp;
    node[u].pos=pos++;
    fp[node[u].pos]=u;
    if(node[u].son==-1) return;
    getpos(node[u].son, sp);
    for(int p=head[u]; p!=-1; p=e[p].next) {
        int v=e[p].to;
        if(v==node[u].son||v==node[u].fa)
            continue;
        getpos(v,v);
    }
}
void build(int id, int l, int r) {
    if(l==r) {
        tr[id] = make_pair(node[fp[l]].v, node[fp[l]].v);
        return;
    }
    int mid=(l+r)/2;
    build(id<<1, l, mid);
    build(id<<1|1, mid+1, r);
    tr[id] = tr[id<<1]*tr[id<<1|1];
}
void update(int id, int l, int r, int pos, int v) {
    if(l==r) {
        tr[id] = make_pair(v,v);
        return; 
    }
    int mid=(l+r)/2;
    if(pos<=mid) update(id<<1, l, mid, pos, v);
    else update(id<<1|1, mid+1, r, pos, v);
    tr[id] = tr[id<<1]*tr[id<<1|1];
}
pii query(int id, int l, int r, int al, int ar) {
    if(al<=l&&r<=ar) return tr[id];
    pii ret=make_pair(-0x3f3f3f3f,0);
    int mid=(l+r)/2;
    if(al<=mid&&ar>=l) ret = ret*query(id<<1, l, mid, al, ar);
    if(al<=r&&ar>mid) ret=ret*query(id<<1|1, mid+1, r, al, ar);
    return ret;
}
pii findans(int u, int v) {
    int f1=node[u].top, f2=node[v].top;
    pii ret=make_pair(-0x3f3f3f3f,0), tmp;
    while(f1!=f2) {
        if(node[f1].dep<node[f2].dep) {
            swap(f1, f2);
            swap(u, v);
        }
        ret = ret*query(1,1,n,node[f1].pos,node[u].pos);
        u=node[f1].fa;
        f1=node[u].top;
    }
    if(node[u].dep>node[v].dep) swap(u,v);
    ret = ret*query(1,1,n,node[u].pos, node[v].pos);
    return ret;
}
int main() {
    scanf("%d", &n);
    init();
    for(int i=0; i<n-1; ++i) {
        int u, v;
        scanf("%d%d", &u, &v);
        e[i<<1]=Edge(v, head[u]);
        head[u] = i<<1;
        e[i<<1|1]=Edge(u, head[v]);
        head[v] = i<<1|1;
    }
    for(int i=1; i<=n; ++i)
        scanf("%d", &node[i].v);
    dfs1(1,0,0);
    getpos(1,1);
    build(1,1,n);
    int q; scanf("%d", &q);
    while(q--) {
        char op[12];
        int u, v;
        scanf("%s%d%d", op, &u, &v);
        if(op[1]=='H')  {
            node[u].v=v;
            update(1,1,n, node[u].pos, v);
        }
        else {
            pair<int,int> ans=findans(u,v);
            if(op[1]=='M') printf("%d\n", ans.first);
            else printf("%d\n", ans.second);
        }
    }
    return 0;
}

//https://vjudge.net/contest/169617#problem/A
#include<cstdio>
#include<cstring>
#include<cmath>
const int MAXN = 128, MAXM=1e4+1e3;
const double INF = 1e10 ;
double p[MAXN][2];
double in[MAXN];
int pre[MAXN], id[MAXN], col[MAXN];
int dbg=0;
struct Edge{
    int u, v;
    double d;
}e[MAXM];
int n, m;
double squ(double x) {
    return 1.0*x*x;
}
void printe() {
    for(int i=1; i<=m; ++i)
        if(e[i].u!=e[i].v)
            printf("%d,%d:%.3f\n",
                    e[i].u, e[i].v, e[i].d);
}
void printin() {
    for(int i=1; i<=n; ++i) {
        printf("%.3f ", in[i]);
    }
    puts("");
}
double ZLA(int root) {
    double ret=0;
    while(true) {
        if(dbg) printe();
        for(int i=1; i<=n; ++i)
            in[i] = INF;
        in[root] = 0;
        for(int i=1; i<=m; ++i) {
            int u=e[i].u, v=e[i].v;
            double d=e[i].d;
            if(u!=v&&in[v]>d) {
                pre[v]=u;
                in[v] = d;
            }
        }
        if(dbg) printin();
        for(int i=1; i<=n; ++i)
            if(in[i]==INF)
                return -1;
        int cnt=0;
        memset(col, -1, sizeof(col));
        memset(id, -1, sizeof(id));
        for(int i=1; i<=n; ++i)
        {
            ret += in[i];
            int v=i;
            while(col[v]!=i&&id[v]==-1&&v!=root) {
                col[v] = i;
                v=pre[v];
            }
            if(v!=root&&id[v]==-1) {
                for(int j=pre[v]; j!=v; j=pre[j])
                    id[j] = cnt+1;
                id[v] = ++cnt;
            }
        }
        if(!cnt) break;
        for(int i=1; i<=n; ++i)
            if(id[i]==-1) id[i] = ++cnt;
        for(int i=1; i<=m; ++i) {
            int v=e[i].v;
            e[i].u=id[e[i].u];
            e[i].v=id[e[i].v];
            if(e[i].u!=e[i].v) e[i].d-=in[v];
        }
        root = id[root];
        n=cnt;
    }
    return ret;
}
int main() {
    while(~scanf("%d%d", &n, &m)) {
        for(int i=1; i<=n; ++i)
            scanf("%lf%lf", p[i], p[i]+1);
        for(int i=1; i<=m; ++i) {
            int u, v;
            scanf("%d%d", &u, &v);
            double d = sqrt(squ(p[u][0]-p[v][0])
                    +squ(p[u][1]-p[v][1]));
            e[i] = (Edge){u,v,d};
        }
        double res = ZLA(1);
        if(res<0) puts("poor snoopy");
        else printf("%.2f\n", res);
    }
    return 0;
}

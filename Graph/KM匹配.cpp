#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
const int MAXN=304;
int slack[MAXN], link[MAXN];
int lx[MAXN], ly[MAXN];
int a[MAXN][MAXN], n;
bool S[MAXN], T[MAXN];
bool dfs(int x) {
    S[x]=true;
    for(int y=0; y<n; ++y) {
        if(T[y]) continue;
        int tmp = lx[x]+ly[y]-a[x][y];
        if(tmp==0) {
            T[y] = true;
            if(link[y]==-1||dfs(link[y])) {
                link[y] = x;
                return true;
            }
        }
        else if(slack[y]>tmp)
            slack[y] = tmp;
    }
    return false;
}
void KM() {
    memset(lx, 0, sizeof(lx));
    memset(ly, 0, sizeof(ly));
    memset(link, -1, sizeof(link));
    for(int x=0; x<n; ++x)
        for(int y=0; y<n; ++y)
            lx[x] = max(lx[x], a[x][y]);
    for(int x=0; x<n; ++x) {
        memset(slack, 0x3f, sizeof(slack));
        while(true) {
            memset(S, false, sizeof(S));
            memset(T, false, sizeof(T));
            if(dfs(x)) break;
            int d=0x3f3f3f3f;
            for(int y=0; y<n; ++y) if(!T[y])
                    d = min(d, slack[y]);
            for(int i=0; i<n; ++i) if(S[i])
                    lx[i] -= d;
            for(int y=0; y<n; ++y) {
                if(T[y]) ly[y] += d;
                else slack[y] -= d;
            }

        }
    }
}
int main() {
    while(~scanf("%d", &n)) {
        for(int i=0; i<n; ++i)
            for(int j=0; j<n; ++j)
                scanf("%d", a[i]+j);
        KM();
        int ans=0;
        for(int i=0; i<n; ++i)
            ans += lx[i]+ly[i];
        printf("%d\n", ans);
    }
    return 0;
}

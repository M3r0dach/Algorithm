#include <cstdio>//don't use this in poj, fzu, zoj
#include <cstring>
#include <queue>
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXL = 1e6+1e3;
const int MAXN = 512;
char s[MAXL], t[MAXN];
int ans[MAXN], tosame[MAXN];
struct Trie{
    int sz, ch[MAXN*MAXN][26];
    int f[MAXN*MAXN], val[MAXN*MAXN];
    int last[MAXN*MAXN];
    void clear(int x) {
        CLR(ch[x]);
        val[x]=f[x]=last[x]=0;
    }
    void init() {
        sz=1;
        clear(0);
    }
    void insert(char *s, int v) {
        int u=0;
        for(char* p=s ; *p; ++p) {
            int idx=*p-'a';
            if(!ch[u][idx]) {
                ch[u][idx]=sz;
                clear(sz++);
            }
            u=ch[u][idx];
        }
        if(val[u]) tosame[v]=val[u];
        else val[u]=v;
    }
    void doc(int u) {
        ++ans[val[u]];
        if(last[u]) doc(last[u]);
    }
    void find(char *s) {
        int u=0;
        for(char* p=s ; *p; ++p) {
            int idx=*p-'a';
            u=ch[u][idx];
            if(val[u]) doc(u);
            else if(last[u]) doc(last[u]);
        }
    }
    void getFail() {
        queue<int> q;
        while (!q.empty()) {
            q.pop();
        }
        f[0]=0;
        rep(i, 0, 26) {
            int u=ch[0][i];
            if(!u) continue;
            q.push(u);
            last[u]=f[u]=0;
        }
        while (!q.empty()) {
            int r=q.front();
            q.pop();
            rep(i, 0, 26) {
                int u=ch[r][i];
                if(!u) {
                    ch[r][i]=ch[f[r]][i];
                    continue;
                }
                q.push(u);
                f[u] = ch[f[r]][i];
                last[u]=val[f[u]]?f[u]:last[f[u]];
            }
        }
    }
}trie;
void solve() {
    int n;
    sf("%d%s", &n, s);
    trie.init();
    CLR(ans);
    CLR(tosame);
    rep(i, 1, n+1) {
        sf("%s", t);
        trie.insert(t, i);
    }
    trie.getFail();
    trie.find(s);
    rep(i, 1, n+1) {
        if(tosame[i]) ans[i]=ans[tosame[i]];
        pf("%d\n", ans[i]);
    }
}
int main() {
    int t=1, ca=0;
    sf("%d", &t);
    while(t--) {
        pf("Case %d:\n", ++ca);
        solve();
    }
    return 0;
}

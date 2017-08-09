#include<cstdio>
#include<cstring>
#include<algorithm>
#include<vector>
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = 2e5+1e3;
struct state {
    int len, pre, ch[26];
    ll cnt;
};
struct SAM {
    int sz, last;
    state st[MAXN];
    state& operator[] (int x) {
        return st[x];
    }
    SAM() {
        sz=1, last=0;
        st[0].len=0, st[0].pre=-1;
    }
    void add(int c) {
        int cur=sz++, p;
        st[cur].len=st[last].len+1;
        for(p=last; p!=-1&&!st[p].ch[c]; p=st[p].pre)
            st[p].ch[c]=cur;
        if(p==-1) st[cur].pre=0;
        else {
            int q=st[p].ch[c];
            if(st[q].len==st[p].len+1)
                st[cur].pre=q;
            else {
                int clone=sz++;
                st[clone]=st[q];
                st[clone].len=st[p].len+1;
                st[cur].pre=st[q].pre=clone;
                for(; p!=-1&&st[p].ch[c]==q; p=st[p].pre)
                    st[p].ch[c]=clone;
            }
        }
        last=cur;
    }
    int find(char *t) {//查询lcs
        int now=0, l=0, ans=0;
        for(char* p=t; *p; ++p) {
            while(now&&!st[now].ch[*p-'a']) {
                now=st[now].pre;
                l=st[now].len;
            }
            if(st[now].ch[*p-'a']) {
                ++l;
                now=st[now].ch[*p-'a'];
            }
            ans=max(l, ans);
        }
        return ans;
    }
} sam;
char s[MAXN];
ll dfs(int u) {
    ll& ret=sam[u].cnt;
    if(ret) return ret;
    ret = 1;
    rep(i, 0, 26)
    if(sam[u].ch[i])
        ret += dfs(sam[u].ch[i]);
    return ret;
}
void printstr(int now, int k) {//输出第K大
    while(--k) {
        rep(i, 0, 26) {
            int v=sam[now].ch[i];
            if(v) {
                if(k<=sam[v].cnt) {
                    putchar(i+'a');
                    now=v;
                    break;
                }
                else k-= sam[v].cnt;
            }
        }
    }
    puts("");
}
int main() {
    sf("%s", s);
    for(char* p=s; *p; ++p)
        sam.add(*p-'a');
    dfs(0);
    int q;
    sf("%d", &q);
    while(q--) {
        int k;
        sf("%d", &k);
        printstr(0, k+1);
    }
    return 0;
}

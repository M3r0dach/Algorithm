//CF 55D
#include<iostream>
#include<cstring>
using namespace std;
typedef long long ll;
const int MAXL = 20;
const int MAXN = 2530;
const int MAXS = 50;
ll dp[MAXL][MAXN][MAXS];
int d[MAXL];
int hashv[MAXN];
void init() {
    int s=2520, cnt=0;
    for(int i=1; i<50; ++i)
        if(s%i==0) {
            hashv[i] = ++cnt;
            hashv[s/i] = 49-cnt;
        }
}
ll gcd(int a, int b) {
    return b?gcd(b,a%b):a;
}
ll dfs(int pos, int r, int lcm, bool lim) {
    if(pos==-1) {
        if(r%lcm==0) return 1;
        return 0;
    }
    int hv = hashv[lcm];
    if(!lim&&dp[pos][r][hv]!=-1)
        return dp[pos][r][hv];
    int up=lim?d[pos]:9;
    ll ret=0;
    for(int i=0; i<=up; ++i) {
        int nr, nlcm;
        if(i) {
            nr = (r*10+i)%2520;
            nlcm = lcm*i/gcd(lcm,i);
        }
        else {
            nr = r*10%2520;
            nlcm = lcm;
        }
        ret += dfs(pos-1, nr, nlcm, lim&&i==up);
    }
    if(!lim) dp[pos][r][hv]=ret;
    return ret;
}
ll getsum(ll x) {
    int pos=0;
    while(x) {
        d[pos++] = x%10;
        x/=10;
    }
    return dfs(pos-1, 0, 1, true);
}
int main() {
    int t;cin >> t;
    init();
    memset(dp, -1, sizeof(dp));
    while(t--) {
        ll l, r;
        cin >> l >> r;
        cout << getsum(r)-getsum(l-1) <<'\n';
    }
    return 0;
}

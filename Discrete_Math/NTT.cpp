#include<cstdio>
#include<cstring>
#include<algorithm>
#include<vector>
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = 1<<19;
const ll MOD = 998244353ll;
ll inv[MAXN];
ll fact[MAXN], ifact[MAXN];
void init() {
    inv[1]=1;
    rep(i, 2, MAXN)
        inv[i] = -(MOD/i)*inv[MOD%i]%MOD;
    fact[0] = ifact[0] = 1;
    rep(i, 1, MAXN) {
        fact[i] = i*fact[i-1]%MOD;
        ifact[i] = inv[i]*ifact[i-1]%MOD;
    }
}
ll c[MAXN], s[MAXN];
ll mypow(ll a, ll x) {
    ll ret=1;
    while(x) {
        if(x&1) ret=ret*a%MOD;
        a=a*a%MOD;
        x >>= 1;
    }
    return ret;
}
int rev(int x, int n) {
    int ret=0;
    while(n>1) {
        ret =ret<<1|(x&1);
        n>>=1;
        x>>=1;
    }
    return ret;
}
void ntt(ll* a, int n, int t) {
    rep(i, 0, n) {
        int rv = rev(i,n);
        if(i<rv) swap(a[i], a[rv]);
    }
    ll g=1;
    if(t==1) g=3;
    else g=inv[3];
    for(int m=2; m<n+1; m<<=1) {
        ll wm= mypow(g, (MOD-1)/m);
        int mid=m>>1;
        for(ll *p=a; p<a+n; p+=m) {
            ll w=1;
            ll *a1=p, *a2=p+mid;
            rep(i, 0, mid) {
                ll t=w*(*a2);
                *a2 = (*a1-t)%MOD;
                *a1 = (*a1+t)%MOD;
                w = w*wm%MOD;
                ++a1, ++a2;
            }
        }
    }
}
void fft(ll *a, ll* b, int n) {
    int tn=MAXN;
    while(tn/4>=n) tn>>=1;
    ntt(a, tn, 1);
    ntt(b, tn, 1);
    rep(i, 0, tn)
        a[i] = a[i]*b[i]%MOD;
    ntt(a, tn, -1);
    rep(i, 0, tn)
        a[i]=a[i]*inv[tn]%MOD;
}
int main() {
    int n;
    init();
    while(~sf("%d", &n)) {
        ++n;
        CLR(s); CLR(c);
        drep(i, n-1, -1) {
            sf("%I64d", c+i);
            c[i] = c[i]*fact[n-1-i]%MOD;
        }
        int m;
        sf("%d", &m);
        ll sum=0, tot=1;
        while(m--) {
            ll x;sf("%I64d", &x);
            sum = sum-x;
        }
        sum %= MOD;
        rep(i, 0, n) {
            s[i] = tot*ifact[i]%MOD;
            tot = tot*sum%MOD;
        }
        fft(c,s,n);
        rep(i, 0, n) {
            ll ans=c[n-1-i]*ifact[i]%MOD;
            if(ans<0) ans+=MOD;
            printf("%lld ", ans);
        }
        puts("");
    }
    return 0;
}

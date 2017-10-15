#include <bits/stdc++.h>//don't use this in poj, fzu, zoj
using namespace std;
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
#define dbg(x) std::cout << #x"=" << x <<'\n';;
typedef long long ll;
const int MAXN = 1e5+1e3;
struct Node{
    static ll p, omega;
    ll a, b;
    Node(ll aa, ll bb) {
        a=aa%p, b=bb%p;
    }
    Node operator *(const Node& rhs) {
        return Node(a*rhs.a%p+b*rhs.b%p*omega%p, a*rhs.b%p+b*rhs.a%p);
    }
};
ll Node::p, Node:: omega;
ll pow_mod(ll a, ll x, ll p) {
    ll ret=1;
    while (x) {
        if(x&1) ret=ret*a%p;
        a=a*a%p;
        x>>=1;
    }
    if(ret<0) ret+=p;
    return ret;
}
Node pow_mod(Node a, ll x) {
    Node ret(1,0);
    while(x) {
        if(x&1) ret = ret*a;
        a=a*a;
        x>>=1;
    }
    return ret;
}
int lgd(ll a, ll p) {
    if(pow_mod(a, (p-1)>>1, p)==1)
        return 1;
    return -1;
}
ll modsqrt(ll a,ll p){
    a%=p;
    if(p==2||a==0) return a%p;
    if(lgd(a, p)==-1) return -1;
    ll x;
    //if(p%4==3) x=pow_mod(a,(p+1)>>2,p);
    //else{
        ll a_0;
        for(a_0=1; a_0<p; ++a_0) {
            if(lgd(a_0*a_0-a, p)==-1)
                break;
        }
        Node::p=p;
        Node::omega=(a_0*a_0-a)%p;
        Node tmp=Node(a_0, 1);
        tmp = pow_mod(tmp, (p+1)>>1);
        x=tmp.a;
        dbg(tmp.b);
    //}
    if(x<0) x+=p;
    if(x*2>p) x=p-x;
    return x;
}
void solve() {
    int a, n;
    sf("%d%d", &a, &n);
    ll res=modsqrt(a,n);
    if(res==-1) puts("No root");
    else if(res<n-res) printf("%lld %lld\n", res, n-res);
    else printf("%lld\n", res);
}
int main() {
    int t=1, ca=0;
    sf("%d", &t);
    while(t--) {
        solve();
    }
    return 0;
}

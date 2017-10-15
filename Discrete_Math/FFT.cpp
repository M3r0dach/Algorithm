#include <bits/stdc++.h>//don't use this in poj, fzu, zoj
#define rep(a,b,c) for(int (a)=(b); (a)<(c); ++(a))
#define drep(a,b,c) for(int (a)=(b); (a)>(c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = (1<<18)+10;
const double pi = acos(-1.0);

struct Complex{
    double r, i;
    Complex(double _r=0, double _i=0)
        :r(_r), i(_i){}
    Complex operator + (const Complex& b) {
        return Complex(r+b.r, i+b.i);
    }
    Complex operator - (const Complex& b) {
        return Complex(r-b.r, i-b.i);
    }
    Complex operator * (const Complex& b) {
        return Complex(r*b.r-i*b.i, r*b.i+i*b.r);
    }
    Complex operator / (double b) {
        return Complex(r/b, i/b);
    }
    Complex conj() {
        return Complex(r, -i);
    }
};
char s[MAXN], t[MAXN];
int rev[MAXN];
Complex w[MAXN];
int geth(int x) {
    int ret=0;
    while ((1<<ret) < x)
        ++ret;
    return ret;
}
void fft_prepare(int n, int h) {
    rev[0]=0;
    rep(i, 1, n) {
        rev[i]=rev[i>>1]>>1|((i&1)<<(h-1));
    }
    rep(i, 0, n) {
        w[i]=Complex(cos(2*pi*i/n), sin(2*pi*i/n));
        //w[i+(n>>1)]=Complex(-w[i].r, -w[i].i);
    }
    w[n]=Complex(1,0);
}
void dft(Complex *a, int n, int h, int on) {
    rep(i, 0, n) {
        if(i<rev[i]) swap(a[i], a[rev[i]]);
    }
    for(int s=1 ; s<h+1; ++s) {
        int m=1<<s;
        int mid=m>>1;
        int stp=n/m;
        for(Complex *p=a; p<a+n; p+=m) {
            for(int i=0; i<mid; ++i) {
                Complex t=p[mid+i];
                if (on==1) t=t*w[i*stp];
                else t=t*w[n-i*stp];
                p[mid+i] = p[i]-t;
                p[i] = p[i]+t;
            }
        }
    }
    if (on==-1) {
        rep(i, 0, n) {
            a[i]=a[i]/n;
        }
    }
}
void solve() {
    static Complex a[MAXN], b[MAXN];
    static int out[MAXN];
    int n=strlen(s);
    int m=strlen(t);
    int h=geth(n+m);
    int tn=1<<h;
    fft_prepare(tn,h);
    rep(i, 0, tn) {
        b[i]=Complex();
    }
    rep(i, 0, n) {
        b[i].r = s[n-i-1]-'0';
    }
    rep(i, 0, m) {
        b[i].i = t[m-i-1]-'0';
    }
    dft(b, tn, h, 1);
    rep(i, 0, tn) {
        int j=(tn-i)&(tn-1);
        a[i] = (b[i]*b[i]-b[j].conj()*b[j].conj())*Complex(0, -0.25);
    }
    dft(a, tn, h, -1);
    int len=0, x=0;
    CLR(out);
    rep(i, 0, tn) {
        x=int(x+a[i].r+0.5);
        out[i]=x%10;
        x/=10;
        if(out[i]) len=i+1;
    }
    while (x) {
        out[len++]=x%10;
        x/=10;
    }
    while (len>1&&out[len-1]==0) {
        --len;
    }
    if (!len) len=1;
    drep(i, len-1, -1)
        printf("%d", out[i]);
    puts("");
}
int main() {
    while(~sf("%s%s", s, t)) {

        solve();
    }
    return 0;
}

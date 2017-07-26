//HDU 4609
#include<cstdio>
#include<cstring>
#include<complex>
#include<algorithm>
#define FIN freopen("input.txt", "r", stdin);
#define FOUT freopen("out.txt", "w", stdout);
using namespace std;
const int MAXN = 1e5+1e3;
const int MAXB = 262144;
const double PI = 3.14159265358979;
typedef long long LL;
complex <double> w[MAXB];
complex <double> fft_pow[MAXB];
void DFT(const complex<double>* a,
         const int& n, const int& step,
         complex<double> *out,
         const int& type) {
    if(n==1) {
        out[0] = a[0];
        return;
    }
    int m = n>>1;
    DFT(a, m, step<<1, out, type);
    DFT(a+step, m, step<<1, out+m, type);
    for(int i=0; i<m; ++i) {
        complex<double> even = out[i];
        complex<double> odd = out[i+m];
        if(type==1)
            odd *= fft_pow[i*step];
        else odd /= fft_pow[i*step];
        out[i] = even+odd;
        out[i+m] = even-odd;
    }
}
void FFT(LL * a, int& n) {
    static complex<double> ta[MAXB];
    static complex<double> tb[MAXB];
    int tn, s;
    for(tn=1; tn<=n+n; tn<<=1);
    for(int i=0; i<tn; ++i)
       if(i<=n) ta[i] = a[i];
       else ta[i] = 0;
    n = tn; s= MAXB/tn;
    for(int i=0; i<n; ++i)
        fft_pow[i] = w[i*s];
    DFT(ta, n, 1, tb, 1);
    for(int i=0; i<n; ++i)
        tb[i] *= tb[i];
    DFT(tb, n, 1, ta, -1);
    for(int i=0; i<n; ++i)
        a[i] = LL(ta[i].real()/n+0.5);
}
int main() {
    int T, n, len;
    static int a[MAXN];
    static LL cnt, total, b[MAXB];
    //FIN;//FOUT;
    for(int i=0; i<MAXB; ++i)
        w[i] = complex<double>(cos(2*PI*i/MAXB), sin(2*PI*i/MAXB));
    scanf("%d", &T);
    while(T--) {
        scanf("%d", &n);
        memset(b, 0, sizeof(b));
        len = 0;
        for(int i=0; i<n; ++i) {
            scanf("%d", a+i);
            len = len>a[i]?len:a[i];
            ++b[a[i]];
        }
        sort(a, a+n);
        FFT(b, len);
        for(int i=0; i<n; ++i)
            --b[a[i]+a[i]];
        for(int i=1; i<len; ++i)
            b[i] = b[i]/2 + b[i-1];
        cnt = 0;
        for(int i=0; i<n; ++i) {
            cnt += b[len-1]-b[a[i]];
            cnt -= (LL)(n-i-1)*i;
            cnt -= n-1;
            cnt -= (LL)(n-i-1)*(n-i-2)/2;
        }
        total = (LL)n*(n-1)/2*(n-2)/3;
        long double res = (long double)(cnt)/total;
        printf("%.7lf\n", (double)res);
    }
    return 0;
}

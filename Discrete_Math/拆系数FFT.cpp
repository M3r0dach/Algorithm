#include <bits/stdc++.h> //don't use this in poj, fzu, zoj
#define rep(a, b, c) for (int(a) = (b); (a) < (c); ++(a))
#define drep(a, b, c) for (int(a) = (b); (a) > (c); --(a))
#define CLR(x) memset(x, 0, sizeof(x))
#define sf scanf
#define pf printf
typedef long long ll;
using namespace std;
const int MAXN = (1 << 17) + 10;
const int MOD = 1e9 + 7;
const double pi = acos(-1.0);
struct Complex {
  double r, i;
  Complex(double _r = 0, double _i = 0) : r(_r), i(_i) {}
  Complex operator+(const Complex &b) { return Complex(r + b.r, i + b.i); }
  Complex operator-(const Complex &b) { return Complex(r - b.r, i - b.i); }
  Complex operator*(const Complex &b) {
    return Complex(r * b.r - i * b.i, r * b.i + i * b.r);
  }
  Complex conj() {
      return Complex(r, -i);
  }
};
ll a[MAXN], b[MAXN], inv[MAXN], c[MAXN];
int geth(int x) {
  int ret = 0, tot = 1;
  while (tot < x) {
    tot <<= 1;
    ++ret;
  }
  return ret;
}
Complex w[MAXN];
int rev[MAXN];
void fft_prepare(int n, int h) {
  rep(i, 0, n >> 1) {
    w[i] = Complex(cos(2 * pi * i / n), sin(2 * pi * i / n));
    w[i + (n >> 1)] = Complex(-w[i].r, -w[i].i);
  }
  w[n] = w[0];
  rev[0] = 0;
  rep(i, 0, n) {
    rev[i] = rev[i >> 1] >> 1;
    if (i & 1)
      rev[i] |= 1 << (h - 1);
  }
}
void dft(Complex *a, int n, int h, int on) {
  rep(i, 0, n) {
    if (i < rev[i])
      swap(a[i], a[rev[i]]);
  }
  for (int m = 2; m < n + 1; m *= 2) {
    int stp = n / m, mid = m >> 1;
    Complex t;
    for (Complex *p = a; p < a + n; p += m) {
      rep(i, 0, mid) {
        if (on > 0)
          t = w[i * stp] * p[mid + i];
        else
          t = w[(n - i * stp) & (n - 1)] * p[mid + i];
        p[mid + i] = p[i] - t;
        p[i] = p[i] + t;
      }
    }
  }
  if (on < 0) {
    rep(i, 0, n) a[i] = Complex(a[i].r / n, a[i].i / n);
  }
}
void conv(ll *a, ll *b, ll* z, int n) {
  int h = geth(n * 2);
  int tn = 1 << h;
  fft_prepare(tn, h);
  static Complex wa[MAXN], wb[MAXN];
  static Complex wc[MAXN], wd[MAXN];
  rep(i, 0, tn) {
    if (i < n) {
      wa[i] = Complex(a[i] >> 15, a[i] & 32767);
      wb[i] = Complex(b[i] >> 15, b[i] & 32767);
    } else
      wb[i] = wa[i] = Complex();
  }
  dft(wa, tn, h, 1);
  dft(wb, tn, h, 1);
  rep(i, 0, tn) {
      int j=(tn-i)&(tn-1);
      Complex dfta = (wa[i]+wa[j].conj())*Complex(0.5,0);
      Complex dftx = (wa[i]-wa[j].conj())*Complex(0, -0.5);
      Complex dftb = (wb[i]+wb[j].conj())*Complex(0.5,0);
      Complex dfty = (wb[i]-wb[j].conj())*Complex(0, -0.5);
      wc[i] = dfta*dftb+dfta*dfty*Complex(0,1);
      wd[i] = dftb*dftx+dftx*dfty*Complex(0,1);
  }
  dft(wc, tn, h, -1);
  dft(wd, tn, h, -1);
  rep(i, 0, tn) {
      ll a=ll(wc[i].r+0.5)%MOD;
      ll b=ll(wc[i].i+0.5)%MOD;
      ll c=ll(wd[i].r+0.5)%MOD;
      ll d=ll(wd[i].i+0.5)%MOD;
      z[i] = ((a<<30)+((b+c)<<15)+d)%MOD;
      if(z[i]<0) z[i] += MOD;
  }
}
void solve() {
  int n, k;
  sf("%d%d", &n, &k);
  rep(i, 0, n) { sf("%I64d", a + i); }
  b[0] = inv[1] = 1;
  rep(i, 2, n) { inv[i] = -(MOD / i) * inv[MOD % i] % MOD; }
  rep(i, 1, n) {
      b[i] = b[i - 1] * (k - 1 + i) % MOD * inv[i] % MOD;
      if(b[i]<0) b[i] += MOD;
  }
  conv(a, b, c, n);
  rep(i, 0, n) {
     printf("%lld\n", c[i]);
  }
}
int main() {
  int t = 1, ca = 0;

  while (t--) {

    solve();
  }
  return 0;
}

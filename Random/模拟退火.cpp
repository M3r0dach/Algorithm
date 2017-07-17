//POJ 1379 run away
#include<cstdio>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#define fi first
#define se second
using namespace std;
const int MAXN=1024;
const int MAXP=10;
const double eps = 1e-4;
const double pi = acos(-1.0);
int X, Y, M;
pair<int,int> p[MAXN];
struct Ans{
    double x, y, dist;
}ans[MAXP];

double sqr(double x) {
    return x*x;
}
double mindis(double x, double y) {
    double ret=1e10;
    for(int i=0; i<M; ++i)
        ret = min(ret, sqr(p[i].fi-x)+sqr(p[i].se-y));
    return ret;
}
double myrand() {
    return (rand()+1.0)/RAND_MAX;
}
void anneal() {
    double d=max(X,Y)/sqrt(1.0*M), r=0.98;
    for(int i=0; i<MAXP; ++i) {
        ans[i].x = myrand()*X;
        ans[i].y = myrand()*Y;
        ans[i].dist = mindis(ans[i].x, ans[i].y);
    }
    while(d>eps) {
        for(int i=0; i<MAXP; ++i) {
            for(int j=0; j<12; ++j) {
                double theta = j*pi/6;
                double nx=ans[i].x+d*cos(theta);
                double ny=ans[i].y+d*sin(theta);
                if(nx>=0&&nx<=X&&ny>=0&&ny<=Y) {
                double dE = mindis(nx,ny)-ans[i].dist;
                if(dE>=0||myrand()<exp(dE/d))
                    ans[i] = (Ans){nx,ny,dE+ans[i].dist};
                }
            }
        }
        d*=r;
    }
    int sel=0;
    for(int i=1; i<MAXP; ++i)
        if(ans[i].dist>ans[sel].dist)
            sel = i;
    printf("The safest point is (%.1f, %.1f).\n", ans[sel].x, ans[sel].y);
}
int main() {
    srand(time(NULL));
    int t;
    scanf("%d", &t);
    while(t--) {
        scanf("%d%d%d", &X, &Y, &M);
        for(int i=0; i<M; ++i)
            scanf("%d%d", &p[i].fi, &p[i].se);
        anneal();
    }
    return 0;
}

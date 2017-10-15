#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
using namespace std;
const double eps=1e-8;
const double PI=acos(-1.0);
struct Point
{
    double x,y;
    Point() {}
    Point(double a,double b)
    {
        x=a,y=b;
    }
    Point operator + (const Point &b)const
    {
        return Point(x+b.x,y+b.y);
    }
    Point operator - (const Point &b)const
    {
        return Point(x-b.x,y-b.y);
    }
    double operator ^(const Point &b)const
    {
        return x*b.y-y*b.x;
    }
    double operator *(const Point &b)const
    {
        return x*b.x+y*b.y;
    }
    bool operator == (const Point &b)const
    {
        return x==b.x&&y==b.y;
    }
};
struct Line
{
    Point s,e;
    double ang;
    Line(){}
    Line(Point a,Point b)
    {
        s=a,e=b;
        ang=atan2(b.y-a.y,b.x-a.x);
    }
    Point operator &(const Line &b)const
    {
        Point res=s;
        double t=((s-b.s)^(b.s-b.e))/((s-e)^(b.s-b.e));
        res.x+=(e.x-s.x)*t;
        res.y+=(e.y-s.y)*t;
        return res;
    }
};
double get_dis(Point a,Point b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
bool Point_on_segment(Point q,Point s,Point e)
{
    if(fabs((q-s)^(e-s))<eps)
        if(min(s.x,e.x)<=q.x&&q.x<=max(s.x,e.x)&&min(s.y,e.y)<=q.y&&q.y<=max(s.y,e.y))
            return 1;
    return 0;
}
bool Segment_cross(Line a,Line b)
{
    if(max(a.e.x,a.s.x)>=min(b.s.x,b.e.x)&&max(b.s.x,b.e.x)>=min(a.s.x,a.e.x)&&max(a.s.y,a.e.y)>=min(b.s.y,b.e.y)&&max(b.s.y,b.e.y)>=min(a.s.y,a.e.y))
        if(((a.s-b.s)^(b.e-b.s))*((b.e-b.s)^(a.e-b.s))>-eps&&((b.s-a.s)^(a.e-a.s))*((a.e-a.s)^(b.e-a.s))>-eps)
            return 1;
    return 0;
}
bool Segment_Line_cross(Line a,Line b)
{
    return ((a.s-b.s)^(b.e-b.s))*((b.e-b.s)^(a.e-b.s))>-eps;
}
Point Point_to_Line(Point P,Line L)
{
    Point result;
    double t = ((P-L.s)*(L.e-L.s))/((L.e-L.s)*(L.e-L.s));
    result.x = L.s.x + (L.e.x-L.s.x)*t;
    result.y = L.s.y + (L.e.y-L.s.y)*t;
    return result;
}
Point Point_to_segment(Point P,Line L)
{
    Point result;
    double t = ((P-L.s)*(L.e-L.s))/((L.e-L.s)*(L.e-L.s));
    if(t>-eps&&t-1<-eps)
    {
        result.x = L.s.x + (L.e.x-L.s.x)*t;
        result.y = L.s.y + (L.e.y-L.s.y)*t;
    }
    else
        result=get_dis(P,L.s)-get_dis(P,L.e)<-eps?L.s:L.e;
    return result;
}
double Cal_area(Point p[],int n)
{
    double res=0;
    for (int i=1;i<n-1;i++)
        res+=(p[i]-p[0])^(p[i+1]-p[0]);
    return res/2.0;
}
Point Rotate(Point vec,double rad)//Anti-clockwise
{
    return Point(vec.x*cos(rad)-vec.y*sin(rad),vec.x*sin(rad)+vec.y*cos(rad));
}
Point Normal_vector(Point vec)
{
    double l=sqrt(vec*vec);
    /*if(l>-eps&&l<eps)
        puts("Zero vector");*/
    return Point(-vec.y/l,vec.x/l);
}
int Point_in_poly(Point p,Point poly[],int n)
{
    int cnt=0;
    Line ray,side;
    ray.s=p,ray.e.y=p.y,ray.e.x=-INF;
    for (int i=0;i<n;i++)
    {
        side.s=poly[i],side.e=poly[(i+1)%n];
        if(Point_on_segment(p,side.s,side.e))
            return 0;
        if(fabs(side.s.y-side.e.y)<eps)
            continue;
        if(Point_on_segment(side.s,ray.s,ray.e)&&side.s.y-side.e.y>eps)
            cnt++;
        else if(Point_on_segment(side.e,ray.s,ray.e)&&side.e.y-side.s.y>eps)
            cnt++;
        else if(Segment_cross(ray,side))
            cnt++;
    }
    return cnt&1?1:-1;
}
int Segment_in_poly(Line a,Point poly[],int n)
{
    if(Point_in_poly(a.s,poly,n)==0||Point_in_poly(a.e,poly,n))
        return 0;
    for (int i=0;i<n;i++)
    {
        if(Segment_Line_cross(a,Line(poly[i],poly[(i+1)%n])))
            return 0;
        if(Point_on_segment(poly[i],a.s,a.e)&&Segment_Line_cross(Line(poly[(i+1)%n],poly[(i+n-1)%n]),a))
            return 0;
    }
    return 1;
}
bool Is_convex(Point poly[],int n)
{
    bool s[3]={0,0,0};
    for (int i=0;i<n;i++)
    {
        int pos;
        if(fabs((poly[(i+1)%n]-poly[i])^(poly[(i+2)%n]-poly[i]))<eps)
            pos=0;
        else
            pos=((poly[(i+1)%n]-poly[i])^(poly[(i+2)%n]-poly[i]))>eps?1:-1;
        s[pos+1]=1;
        if(s[0]&&s[2])
            return 0;
    }
    return 1;
}
bool Point_cmp(Point a,Point b)
{
    if(a.x==b.x)
        return a.y<b.y;
    return a.x<b.x;
}
int Convex_hull(Point p[],int n,Point ans[])
{
    sort(p,p+n,Point_cmp);
    n=unique(p,p+n)-p;
    int m=0;
    for (int i=0;i<n;i++)
    {
        while(m>1&&((ans[m-1]-ans[m-2])^(p[i]-ans[m-2]))<=0)
            m--;
        ans[m++]=p[i];
    }
    int k=m;
    for (int i=n-2;i>=0;i--)
    {
        while(m>k&&((ans[m-1]-ans[m-2])^(p[i]-ans[m-2]))<=0)
            m--;
        ans[m++]=p[i];
    }
    if(n>1)
        m--;
    return m;
}
double Rotating_Calipers(Point poly[],int n)//Anti-clockwise
{
    int j=1;
    double ans=0;
    poly[n]=poly[0];
    for (int i=0;i<n;i++)
    {
        while(((poly[i]-poly[i+1])^(poly[i]-poly[j]))-((poly[i]-poly[i+1])^(poly[i]-poly[j+1]))<-eps)
            j=(j+1)%n;
        ans=max(ans,max(get_dis(poly[i],poly[j]),get_dis(poly[i+1],poly[j])));
    }
    return ans;
}
bool HPI_cmp(Line a,Line b)
{
    if(fabs(a.ang-b.ang)>eps)
        return a.ang<b.ang;
    return ((a.s-b.s)^(b.e-b.s))<-eps;//保留逆时针方向那个
}
int HPI(Line line[],int n,Point ans[])
{
    Line Q[MAXN];
    int tot=1,cnt=0;
    sort(line,line+n,HPI_cmp);
    for (int i=1;i<n;i++)
        if(fabs(line[i].ang-line[i-1].ang)>eps)
            line[tot++]=line[i];
    int head=0,tail=1;
    Q[0]=line[0];
    Q[1]=line[1];
    for(int i = 2; i < tot; i++)
    {
        /*if(fabs((Q[tail].e-Q[tail].s)^(Q[tail-1].e-Q[tail-1].s))<eps||fabs((Q[head].e-Q[head].s)^(Q[head+1].e-Q[head+1].s))<eps)
            return 0;*/
        while(head < tail && (((Q[tail]&Q[tail-1])-line[i].s)^(line[i].e-line[i].s)) > eps)
            tail--;
        while(head < tail && (((Q[head]&Q[head+1])-line[i].s)^(line[i].e-line[i].s)) > eps)
            head++;
        Q[++tail] = line[i];
    }
    while(head < tail && (((Q[tail]&Q[tail-1])-Q[head].s)^(Q[head].e-Q[head].s)) > eps)
        tail--;
    while(head < tail && (((Q[head]&Q[head-1])-Q[tail].s)^(Q[tail].e-Q[tail].e)) > eps)
        head++;
    if(tail<=head+1)
        return 0;
    for (int i=head;i<tail;i++)
        ans[cnt++]=Q[i]&Q[i+1];
    if(head<tail-1)
        ans[cnt++]=Q[head]&Q[tail];
    return cnt;
}

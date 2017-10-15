#include<stdio.h>
#include<algorithm>
using namespace std;
#define N 100100
#define rep(a,b,c) for(int a=b; a<c; ++a)
int data[N];

struct node
{
    int v[N];
    int num[N];
} td[31];
void build(int l,int r,int dep)
{
    if(l>=r)return;
    int i,mid=(l+r)>>1,midd=data[mid],ant=mid-l+1;
    // ant保存有多少和sorted[mid]一样大的数进入左孩子
    int ln=l-1,rn=mid;
    for(i=l; i<=r; i++)
        if(td[dep].v[i]<midd)
            ant--;
    for(i=l; i<=r; i++)
    {
        if(i==l)
            td[dep].num[i]=0;
        else
            td[dep].num[i]=td[dep].num[i-1];
        if(td[dep].v[i]<midd||(td[dep].v[i]==midd&&ant))
        {
            if(td[dep].v[i]==midd)--ant;
            td[dep+1].v[++ln]=td[dep].v[i];
            ++td[dep].num[i];
        }
        else td[dep+1].v[++rn]=td[dep].v[i];
    }
    build(l,mid,dep+1);
    build(mid+1,r,dep+1);
}
int query(int a,int b,int k,int l,int r,int dep)
{
    int mid=(l+r)>>1;
    if(a==b)return td[dep].v[a];
    int lx,rx=td[dep].num[b],sub,sl,sr;
    lx=(a<=l)?0:td[dep].num[a-1];
    sub=rx-lx;
    if(sub>=k)
        return query(l+lx,l+rx-1,k,l,mid,dep+1);
    else
    {
        sl=a-l-lx;
        sr=b-l-rx;
        return query(mid+1+sl,mid+1+sr,k-sub,mid+1,r,dep+1);
    }
}
int n;
void dbg() {
    rep(i, 0, 6) {
        printf("%d\n", i);
        rep(j, 1, n+1)
            printf("%d%c", td[i].num[j], j==n?'\n':' ');
    }
    rep(i, 0, 6) {
        printf("%d\n", i);
        rep(j, 1, n+1)
            printf("%d%c", td[i].v[j], j==n?'\n':' ');
    }
}
int main()
{
    int i,cas=1,m,l,r,k;
  //  scanf("%d",&cas);
    while(cas--)
    {
        scanf("%d%d",&n,&m);
        for(i=1; i<=n; i++)
        {
            scanf("%d",&td[1].v[i]);
            data[i]=td[1].v[i];
        }
        sort(data+1,data+n+1);
        build(1,n,1);
        dbg();
        m=0;
        while(m--)
        {
            scanf("%d%d%d",&l,&r,&k);
            printf("%d\n",query(l,r,k,1,n,1));
        }
    }
    return 0;
}

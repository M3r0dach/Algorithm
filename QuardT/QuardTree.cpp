//HDU-6183
#include <bits/stdc++.h>
using namespace std;
const int maxn = 1<<15;
typedef long long ll;
typedef double db;
#define CLR(x) memset(x, 0, sizeof(x))
const ll mod = 1e9 + 7;
const int Cap = 64;
const int level = 10;
struct Point{
    int x, y, c;
    Point(){}
    Point(int x, int y, int c=0)
    :x(x), y(y), c(c){}
};

struct Rect{
    int x, y, L;
    Rect(int x=0, int y=0, int L=0)
    :x(x),y(y),L(L){}
    int midx() {
        return x+L/2;
    }
    int midy() {
        return y+L/2;
    }
    bool in(Point p) {
        return x<=p.x&&p.x<x+L&&y<=p.y&&p.y<y+L;
    }
};
struct Area{
    int x1, y1, y2;
    Area(int x1=0, int y1=0, int y2=0)
    :x1(x1), y1(y1), y2(y2){}
    bool in(const Point& p) {
        return 1<=p.x&&p.x<=x1&&y1<=p.y&&p.y<=y2;
    }
    bool insec(const Rect& b) {
        int l=max(1, b.x), r=min(x1, b.x+b.L-1);
        int s=max(y1, b.y), h=min(y2, b.y+b.L-1);
        return l<=r&&s<=h;
    }
    bool have(const Rect& b) {
        return 1<=b.x&&b.x+b.L-1<=x1&&y1<=b.y&&b.y+b.L-1<=y2;
    }
};
struct Node{
    Rect bound;
    vector<Point> P;
    bitset<64> st;
    int ch[4], sz;
    void split();
    void insert(Point& p);
    void find(Area& A);
}node[maxn];
int tot;

int newnode(const Rect& rt) {
    node[tot].bound=rt;
    node[tot].P.clear();
    node[tot].sz=0;
    node[tot].st.reset();
    CLR(node[tot].ch);
    return tot++;
}
void Node::split() {
    ch[0]=newnode(Rect(bound.x, bound.y, bound.L/2));
    ch[1]=newnode(Rect(bound.midx(), bound.y, bound.L/2));
    ch[2]=newnode(Rect(bound.midx(), bound.midy(), bound.L/2));
    ch[3]=newnode(Rect(bound.x, bound.midy(), bound.L/2));
    for(int i=0; i<P.size(); ++i)
        for(int j=0; j<4; ++j)
            node[ch[j]].insert(P[i]);
    P.clear();
}
void Node::insert(Point &p) {
    if(!bound.in(p)) return;
    ++sz;
    st.set(p.c);
    if(!ch[0]) {
        P.push_back(p);
        if(P.size()==Cap&&bound.L>1)
            split();
        return;
    }
    for(int i=0; i<4; ++i)
        node[ch[i]].insert(p);
}
bitset<64> ans;
void Node::find(Area& A) {
    if(!A.insec(bound)) return;
    if(A.have(bound)) ans |= st;
    else if(!ch[0]) {
        for(int i=0; i<P.size(); ++i)
            if(A.in(P[i])) ans.set(P[i].c);
    }
    else {
        for(int i=0; i<4; ++i)
            node[ch[i]].find(A);
    }
}
void init() {
    tot=1;
    newnode(Rect(1,1,1<<20));
}
void Add(int x, int y, int c) {
    Point p=Point(x,y,c);
    node[1].insert(p);
}
int Ask(int x1, int y1, int y2) {
    Area A(x1, y1, y2);
    ans.reset();
    node[1].find(A);
    return ans.count();
}

int main(){
    int cmd;
    while(~scanf("%d", &cmd)) {
        if(cmd==0) init();
        else if(cmd==3) break;
        else {
            int a, b, c;
            scanf("%d%d%d", &a, &b, &c);
            if(cmd==1) Add(a,b,c);
            else printf("%d\n", Ask(a, b, c));
        }
    }
    return 0;
}

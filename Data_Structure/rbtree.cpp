#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<ctime>
#define RED false
#define BLACK true
const int MAXN = 1e5+1e3;
struct Node{
    int key;
    bool color;
    Node *ch[2], *p;
    int cmp(Node* b) {
        return b->key > key;
    }
};
struct RBTree{
    Node pool[MAXN];
    int sz;
    Node *root, *nil;
    RBTree() {
        clear();
    }
    void clear() {
        sz = 0;
        nil = pool;
        nil->color = BLACK;
        root = nil;
        root->p = nil;
        sz = 0;
    }
    Node* newNode(int v) {
        Node *ret = &pool[++sz];
        ret->key = v;
        ret->ch[0] = ret->ch[1] = nil;
        ret->color = RED;
        return ret;
    }
    void rotate(Node* x, int d) { // 0 left, 1 right
        Node *y = x->ch[d^1];
        if(y->ch[d]!=nil) {
            y->ch[d]->p = x;
        }
        x->ch[d^1] = y->ch[d];
        y->p = x->p;
        if(y->p==nil) {
            root = y;
        } else {
            y->p->ch[y->p->cmp(y)] = y;
        }
        x->p = y;
        y->ch[d] = x;
    }
    void fixinsert(Node* x) {
        while(x->p->color==RED) {
            Node* pp = x->p->p;
            Node* uncle = pp->ch[pp->cmp(x)^1];
            if(uncle->color==RED) {
                pp->color=RED;
                uncle->color=BLACK;
                x->p->color=BLACK;
                x = pp;
            } else {
                if(x->p->cmp(x)!=pp->cmp(x->p)) {
                    int d = x->p->cmp(x);
                    x = x->p;
                    rotate(x, d^1);
                }
                x->p->color = BLACK;
                pp->color = RED;
                rotate(pp,pp->cmp(x->p)^1);
            }
        }
        root->color = BLACK;
    }
    void insert(int k) {
        Node *z = newNode(k);
        Node *x = root, *y=nil;
        while(x!=nil) {
            y = x;
            x = x->ch[x->cmp(z)];
        }
        z->p = y;
        if(y==nil) {
            root=z;
            return;
        } else {
            y->ch[y->cmp(z)] = z;
        }
        fixinsert(x);
    }
    Node* find(int k) {
        Node *x = root;
        while(x!=nil&&x->key!=k) {
            if(k>x->key) x=x->ch[1];
            else x=x->ch[0];
        }
        return x;
    }
    void transpant(Node* u, Node* v) {
        if(u->p==nil) {
            root = v;
        } else {
            u->p->ch[u->p->cmp(u)] = v;
        }
        v-> p = u->p; //  nil v also need
    }
    Node* min(Node* x) {
        if(x==nil) return x;
        while(x->ch[0]!=nil) {
            x=x->ch[0];
        }
        return x;
    }
    void fixdel(Node* x) {
        while(x!=root&&x->color==BLACK) {
            int d = x->p->ch[1]==x;
            Node* w= x->p->ch[d^1];
            if(w->color==RED) {
                w->color = BLACK;
                x->p->color = RED;
                rotate(x->p, d);
                w=x->p->ch[d^1];
            }
            if(w->ch[0]->color==BLACK&&w->ch[1]->color==BLACK) {
                w->color=RED;
                x = x->p;
            }
            else {
                if(w->ch[d^1]->color==BLACK) {
                    w->ch[d]->color = BLACK;
                    w->color=RED;
                    rotate(w, d^1);
                    w = x->p->ch[d^1];
                }
                x->p->color = BLACK;
                w->color = RED;
                w->ch[d^1]->color = BLACK;
                rotate(x->p, d);
                x = root;
            }
        }
        x->color = BLACK;
    }
    void remove(Node* z) {
        Node *y=z, *x=nil;
        bool origincolor = y->color;
        if(y->ch[0]==nil) {
            x=y->ch[1];
            transpant(y,x);
        } else if(y->ch[1]==nil) {
            x=y->ch[0];
            transpant(y,x);
        } else {
            y = min(z->ch[1]);
            x = y->ch[1];
            if(y->p!=z) {
                transpant(y,x);
                y->ch[1] = z->ch[1];
                y->ch[1]->p = y;
            }
            y->ch[0] = z->ch[0];
            y->ch[0]->p = y;
            y->color = z->color;
            transpant(z,y);
        }
        if(origincolor==BLACK) {
            fixdel(x);
        }
    }
    void remove(int k) {
        Node* x= find(k);
        if(x!=nil) remove(x);
    }
}rbt;
int a[MAXN], n=MAXN-10;
int main() {
    int sz = n;
    for(int i=0; i<sz; i++) {
        a[i] = i;
    }
    srand(time(0));
    while(sz) {
        int k = rand()%sz;
        rbt.insert(a[k]);
        a[k] = a[sz-1];
        sz--;
    }
    for(int i=1; i<=n; i++) {
        int v = rbt.min(rbt.root)->key;
        printf("%d\n", v);
        rbt.remove(v);
    }
    return 0;
}
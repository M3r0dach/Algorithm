#include<cstdio>
#include<cstring>
const int MAXN = 2048;
char s[MAXN];
int l;
int min(int a, int b) {
    return a<b?a:b;
}
void print(int st) {
    for(int i=0; i<l; ++i)
        putchar(s[i+st]);
    puts("");
}
int main() {
    scanf("%s", s);
    l = strlen(s);
    for(int i=0; i<l; ++i)
        s[i+l] = s[i];
    s[l<<1] = 0;
    if(l==1) {
        print(0);
        return 0;
    }
    int i, j;
    for(i=0, j=1; i<l&&j<l; ) {
        int k=0;
        while(k<l&&s[i+k]==s[j+k]) ++k;
        if(k==l) break;
        if(s[i+k]>s[j+k]) i=i+k+1;
        else j=j+k+1;
        if(i==j) j+=1;
    }
    print(min(i,j));
    return 0;
}

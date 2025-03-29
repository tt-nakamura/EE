#include "EE.h"

std::ostream& operator<<(std::ostream& s, const EE& a) {// output a to s
    s << '[' << a.x << ' ' << a.y << ']';
    return s;
}

void conv(EE& b, const ZZ& a) { b.x = a; clear(b.y); }
void conv(EE& b, long a) { b.x = a; clear(b.y); }

EE& EE::operator=(const ZZ& a) { x=a; clear(y); return *this; }
EE& EE::operator=(long a) { x=a; clear(y); return *this; }

EE& operator+=(EE& a, const EE& b) {// a+=b
    a.x += b.x;
    a.y += b.y;
    return a;
}

EE& operator-=(EE& a, const EE& b) {// a-=b
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

EE& operator*=(EE& a, const EE& b) {// a*=b
    mul(a,a,b);
    return a;
}

EE& operator/=(EE& a, const EE& b) {// a/=b
    div(a,a,b);
    return a;
}

EE& operator%=(EE& a, const EE& b) {// a%=b
    rem(a,a,b);
    return a;
}

void set(EE& a) { set(a.x); clear(a.y); }// a=1
void set(EE& a, const ZZ& x, const ZZ& y) { a.x=x; a.y=y; }// a=x+yw
void set(EE& a, long x, long y) { a.x=x; a.y=y; }// a=x+yw
void clear(EE& a) { clear(a.x); clear(a.y); }// a=0

void conj(EE& b, const EE& a) {// b = complex conjugate of a
    sub(b.x, a.x, a.y);
    negate(b.y, a.y);
}

void mirror(EE& b, const EE& a)
// b = reflection of a with respect to line y=x/sqrt(3)
{
    if(&b!=&a) b.x = a.x;
    sub(b.y, a.x, a.y);
}

void norm(ZZ& x, const EE& a) {// x = |a|^2
    sub(x, a.x, a.y);
    sqr(x,x);
    MulAddTo(x, a.x, a.y);
}

void negate(EE& b, const EE& a) {// b=-a
    negate(b.x, a.x);
    negate(b.y, a.y);
}

void rot60(EE& b, const EE& a) {// b = a*(1+w)
    if(&b==&a) {
        swap(b.x, b.y);
        sub(b.x, b.y, b.x);
    }
    else {
        sub(b.x, a.x, a.y);
        b.y = a.x;
    }
}

void rot120(EE& b, const EE& a) {// b = a*w
    if(&b==&a) {
        swap(b.x, b.y);
        b.y -= b.x;
        negate(b.x, b.x);
    }
    else {
        sub(b.y, a.x, a.y);
        negate(b.x, a.y);
    }
}

void rot60(EE& b, const EE& a, long k) {// b = a*(1+w)^k
    if((k%=6)==1) rot60(b,a);
    else if(k==2) rot120(b,a);
    else if(k==3) negate(b,a);
    else if(k==4) { rot60(b,a); negate(b,b); }
    else if(k==5) { rot120(b,a); negate(b,b); }
    else if(&b!=&a) b=a;
}

long hexant(const EE& a)
// return -1 if a==0
// return k=0,1,...,5 if 60k <= arg(a) < 60(k+1)
{
    if(IsZero(a.y)) {
        if(IsZero(a.x)) return -1;
        else if(sign(a.x) > 0) return 0;
        else return 3;
    }
    else if(sign(a.y) > 0) {
        if(a.x > a.y) return 0;
        else if(sign(a.x) > 0) return 1;
        else return 2;
    }
    else if(a.x < a.y) return 3;
    else if(sign(a.x) < 0) return 4;
    else return 5;
}

long FirstHex(EE& b, const EE& a)
// b = a*(1+w)^k for some k=0,1,...,5
//     in first hexant 0 <= arg(b) < 60
// return k
// if a is zero, then b=0 and k=0
{
    long k(hexant(a));
    if(k>0) k=6-k; else k=0;
    rot60(b,a,k);
    return k;
}

void add(EE& c, const EE& a, const EE& b) {// c=a+b
    add(c.x, a.x, b.x);
    add(c.y, a.y, b.y);
}

void sub(EE& c, const EE& a, const EE& b) {// c=a-b
    sub(c.x, a.x, b.x);
    sub(c.y, a.y, b.y);
}

void mul(EE& c, const EE& a, const EE& b) {// c=a*b
    ZZ s,t,u,v;
    mul(s, a.x, b.x);
    mul(t, a.y, b.y);
    sub(u, a.y, a.x);
    sub(v, b.x, b.y);
    sub(c.x, s, t);
    mul(c.y, u, v);
    c.y += s;
}

void sqr(EE& b, const EE& a) {// b=a*a
    ZZ s,t;
    add(s, a.x, a.y);
    sub(t, a.x, a.y);
    mul(b.x, s, t);
    t <<= 1;
    t += a.y;
    mul(b.y, a.y, t);
}

void div(EE& q, const EE& a, const EE& b)
// q = quotient of a/b such that
//     a = bq + r and norm(r) <= (3/4)norm(b)
{
    ZZ n;
    EE c;
    if(IsZero(b.y)) {
        n = b.x;
        c = a;
    }
    else if(IsZero(b.x)) {
        negate(n, b.y);
        rot60(c,a);
    }
    else if(b.x == b.y) {
        negate(n, b.x);
        rot120(c,a);
    }
    else {
        norm(n,b);
        conj(c,b);
        c *= a;
    }
    c.x <<= 1; c.x += n;
    c.y <<= 1; c.y += n;
    n <<= 1;
    div(q.x, c.x, n);
    div(q.y, c.y, n);
}

void rem(EE& r, const EE& a, const EE& b)
// r = remainder of a/b such that
//     a = bq + r and norm(r) <= (3/4)norm(b)
{
    EE q;
    DivRem(q,r,a,b);
}

void DivRem(EE& q, EE& r, const EE& a, const EE& b)
// q,r = quotient and remainder of a/b such that
//       a = bq + r and norm(r) <= (3/4)norm(b)
{
    if(&a==&q || &a==&r) {
        EE c(a);
        DivRem(q,r,c,b);
        return;
    }
    if(&b==&q) {
        EE c(b);
        DivRem(q,r,a,c);
        return;
    }
    div(q,a,b);
    mul(r,b,q);
    sub(r,a,r);
}

long divide(EE& q, const EE& a, const EE& b)
// if a/b is divisible, set q=a/b and return 1
// else return 0 and q is unchanged
{
    ZZ n;
    EE c;
    if(IsZero(b.y)) {
        n = b.x;
        c = a;
    }
    else if(IsZero(b.x)) {
        negate(n, b.y);
        rot60(c,a);
    }
    else if(b.x == b.y) {
        negate(n, b.x);
        rot120(c,a);
    }
    else {
        norm(n,b);
        conj(c,b);
        c *= a;
    }
    if(!divide(c.x, c.x, n) ||
       !divide(q.y, c.y, n)) return 0;
    q.x = c.x;
    return 1;
}

long divide(const EE& a, const EE& b)
// if a/b is divisible, return 1, else return 0
{
    ZZ n;
    EE c;
    if(IsZero(b.y) || b.x == b.y)
        return divide(a.x, b.x) && divide(a.y, b.x);
    else if(IsZero(b.x))
        return divide(a.x, b.y) && divide(a.y, b.y);
    norm(n,b);
    conj(c,b);
    c *= a;
    return divide(c.x, n) && divide(c.y, n);
}

long divide3(EE& q, const EE& a)
// if a is divisible by 1-w, set q=a/(1-w) and return 1
// else return 0 and q is unchanged
{
    ZZ n;
    EE c;
    rot120(c,a);
    LeftShift(n, a.x, 1); c.x += n;
    LeftShift(n, a.y, 1); c.y += n;
    if(!divide(c.x, c.x, 3) ||
       !divide(q.y, c.y, 3)) return 0;
    q.x = c.x;
    return 1;
}

void power(EE& b, const EE& a, long n)
// b = a^n; assume n>=0
{
    if(n==0 || IsOne(a)) { set(b); return; }
    if(&b==&a) { EE c(a); power(b,c,n); return; }
    long m(1<<(NumBits(n)-1));
    b=a;
    for(m>>=1; m; m>>=1) {
        sqr(b,b);
        if(n&m) b*=a;
    }
}

long IsUnit(const EE& a) {// test if |a|==1
    return (IsZero(a.y) && (IsOne(a.x) || a.x==-1)) ||
        ((IsOne(a.y) || a.y==-1) && (IsZero(a.x) || a.x==a.y));
}

long IsAssoc(const EE& a, const EE& b)
// return 1 if a = b*(1+w)^k for some k=0,1,...,5
// else return 0
{
    EE c,d;
    FirstHex(c,a);
    FirstHex(d,b);
    return c==d;
}

void GCD(EE& d, const EE& a, const EE& b)
// d = greatest common divisor of a and b
//     in first hexant 0 < d.y < d.x
// by euclidean algorithm
{
    EE x(a),y(b),r;
    while(!IsZero(y)) {
        rem(r,x,y);
        x=y;
        y=r;
    }
    FirstHex(d,x);
}

void XGCD(EE& d, EE& s, EE& t, const EE& a, const EE& b)
// d = greatest common divisor of a and b
//     in first hexant 0 < d.y < d.x
//     and compute s,t such that d = s*a + t*b
// by extended euclidean algorithm
{
    long k;
    EE x(a),y(b),u,v(1),q,r;
    set(s);
    clear(t);
    while(!IsZero(y)) {
        DivRem(q,r,x,y);
        mul(x,q,u);
        sub(x,s,x);
        s=u;
        u=x;
        mul(x,q,v);
        sub(x,t,x);
        t=v;
        v=x;
        x=y;
        y=r;
    }
    k = FirstHex(d,x);
    rot60(s,s,k);
    rot60(t,t,k);
}

// b = random Eisenstein integer in hexagon R*e^{i\pi n/3} (n=0...5)
void RandomBits(EE& b, long l)// 0 <= R < 2^l
{ RandomBits(b.x, l); RandomBits(b.y, l); rot60(b, b, RandomBnd(3)<<1); }

void RandomLen(EE& b, long l)// 2^{l-1} <= R < 2^l 
{ RandomLen(b.x, l); RandomLen(b.y, l); rot60(b, b, RandomBnd(3)<<1); }

void RandomBnd(EE& b, const ZZ& a)// 0 <= R < a
{ RandomBnd(b.x, a); RandomBnd(b.y, a); rot60(b, b, RandomBnd(3)<<1); }

void PowerMod(EE& b, const EE& a, long n, const EE& m)
// b = a^n mod m; assume n>=0 and |a| < |m|
{
    if(n==0 || IsOne(a)) { set(b); return; }
    if(&b==&a) { EE c(a); PowerMod(b,c,n,m); return; }
    long k(1<<(NumBits(n)-1));
    b=a;
    for(k>>=1; k; k>>=1) {
        sqr(b,b); b%=m;
        if(n&k) { b*=a; b%=m; }
    }
}

void PowerMod(EE& b, const EE& a, const ZZ& n, const EE& m)
// b = a^n mod m; assume n>=0 and norm(a) < norm(m)
{
    if(IsZero(n) || IsOne(a)) { set(b); return; }
    if(&b==&a) { EE c(a); PowerMod(b,c,n,m); return; }
    b=a;
    for(long k=NumBits(n)-2; k>=0; k--) {
        sqr(b,b); b%=m;
        if(bit(n,k)) { b*=a; b%=m; }
    }
}

long ProbPrime(const EE& a, long NTRY)
// return 1 if either
//   |a| is prime and |a|==2 (mod 3) or
//   norm(a) is prime and norm(a)==1 (mod 3) or
//   norm(a) == 3
// else return 0
// NTRY: number of trials of Miller-Rabin test
// primality is tested probabilistically MR method
{
    ZZ b;
    if(IsZero(a.y) || a.x==a.y) abs(b, a.x);
    else if(IsZero(a.x)) abs(b, a.y);
    if(!IsZero(b))
        return b%3 == 2 && ProbPrime(b, NTRY);
    norm(b,a);
    return b==3 || b%3 == 1 && ProbPrime(b, NTRY);
}

void GenPrime(EE& p, long l, long f, long err)
// generate random Eisenstein prime p
// Assume f=1 or 2; Assume l>=2
// if f=1, |p|^2=q and q==1 (mod 3)
// if f=2, p=q+0i and q==2 (mod 3)
// where q is random prime and 2^{l-1} <= q < 2^l
// probability of error is less than 2^-err
// p is primary, i.e., p.x==2 and p.y==0 (mod 3)
{
    if(l<2) Error("l<2 in GenPrime");
    ZZ q;
    if(f==1) {
        do GenPrime(q,l,err);
        while(q%3 == 2);
        if(q==3) { set(p,1,-1); return; }
        FactorPrime(p,q);
        if(RandomBits_long(1)) conj(p,p);
    }
    else if(f==2) {
        do GenPrime(p.x, l, err);
        while(p.x%3 != 2);
        clear(p.y);
    }
    else Error("f<1 or f>2 in GenPrime");
}

long primary(EE& b, const EE& a)
// b = unit * a such that b.x==2 and b.y==0 (mod 3)
// return k=0,1,2 such that a = \pm omega^k * b
// Assume norm(a) != 0 (mod 3)
{
    long x(a.x%3), y(a.y%3);
    if(x==y) {
        rot120(b,a);
        if(y==2) negate(b,b);
        return 2;
    }
    if(x==0) {
        rot60(b,a);
        if(y==2) negate(b,b);
        return 1;
    }
    if(x==1) negate(b,a);
    else if(&b!=&a) b=a;
    return 0;
}

void ResSymb(EE& s, const EE& a, const EE& b)
// s = cubic residue symbol (a/b)_3 = 0,1,w,w^2
// Assume norm(a) < norm(b) and norm(b) != 0 (mod 3)
// Assume b is primary, but may not be prime
// reference: K. Ireland and M. Rosen
//   "A Classical Introduction to Modern Number Theory" section 9.3
{
    long i,j(0),m,n;
    ZZ M,N;
    EE u(a),v(b),w;
    while(!IsZero(u)) {
        DivRem(M, v.x, 3); M++;
        DivRem(N, v.y, 3);
        m = M%3;
        n = AddMod(N%3, m, 3);
        while(divide3(u,u))
            j = SubMod(j,m,3);// supplementary law for 1-w
        for(i=primary(u,u); i; i--)
            j = AddMod(j,n,3);// supplementary law for units
        rem(w,v,u);
        v = u;
        u = w;
    }
    if(!IsUnit(v)) clear(s);
    else if(j==0) set(s);
    else if(j==1) set(s,0,1);
    else set(s,-1,-1);
}
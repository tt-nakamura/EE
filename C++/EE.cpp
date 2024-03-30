// uses NTL
//   http://www.shoup.net/ntl

#include "EE.h"
using namespace NTL;

EE& EE::operator=(const ZZ& a)
{ x=a; clear(y); return *this; }// a+0w

EE& EE::operator=(long a)
{ x=a; clear(y); return *this; }// a+0w

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
void set(EE& a, long x, long y) { a.x=x; a.y=y; }// a = x+yw
void clear(EE& a) { clear(a.x); clear(a.y); }// a=0

void conj(EE& b, const EE& a)
// b = complex congugate of a
{
    sub(b.x, a.x, a.y);
    negate(b.y, a.y);
}

void mirror(EE& b, const EE& a)
// a and b are symmetric wrt line y=x/sqrt(3)
{
    if(&b!=&a) b.x = a.x;
    sub(b.y, a.x, a.y);
}

void norm(ZZ& n, const EE& a) {// n = a*conj(a)
    sub(n, a.x, a.y);
    sqr(n,n);
    MulAddTo(n, a.x, a.y);
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

void rot60e(EE& b, const EE& a, long e) {// b = a*(1+w)^e
    if((e%=6)==1) rot60(b,a);
    else if(e==2) rot120(b,a);
    else if(e==3) negate(b,a);
    else if(e==4) { rot60(b,a); negate(b,b); }
    else if(e==5) { rot120(b,a); negate(b,b); }
    else if(&b!=&a) b=a;
}

long hexant(const EE& a)
// return floor(arg(a)*3/pi) if a!=0
// return -1 if a==0
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
// b = a*(1+w)^e so that 0<=arg(b)<pi/3
// and return e
{
    long e(hexant(a));
    if(e>0) e=6-e; else e=0;
    rot60e(b,a,e);
    return e;
}

void primary(EE& b, const EE& a)
// b = a*unit == x+yw to make x==2 and y==0(mod 3)
// assume norm(a) != 0 (mod 3)
{
    long x(a.x%3), y(a.y%3);
    if(x==y) {
        rot120(b,a);
        if(y==2) negate(b,b);
    }
    else if(x==0) {
        rot60(b,a);
        if(y==2) negate(b,b);
    }
    else if(x==1) negate(b,a);
    else if(&b!=&a) b=a;
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
//   a = bq + r and |r/b|^2 <= 3/4
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
//   a = bq + r and |r/b|^2 <= 3/4
{
    EE q;
    DivRem(q,r,a,b);
}

void DivRem(EE& q, EE& r, const EE& a, const EE& b)
// q,r = quotient and remainder of a/b such that
//   a = bq + r and |r/b|^2 <= 3/4
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
// else return 0 (and q is unchanged)
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
    if(IsZero(b.x))
        return divide(a.x, b.y) && divide(a.y, b.y);
    norm(n,b);
    conj(c,b);
    c *= a;
    return divide(c.x, n) && divide(c.y, n);
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

long IsUnit(const EE& a) {// test if norm(a)==1
    return (IsZero(a.y) && IsOne(a.x) || a.x==-1) ||
           (IsOne(a.y) && IsZero(a.x) || IsOne(a.x)) ||
           (a.y==-1 && IsZero(a.x) || a.x==-1);
}

long IsAssoc(const EE& a, const EE& b)
// return 1 if a = b*w^e for some e (e=0,1,2,3) else 0
{
    EE c,d;
    FirstHex(c,a);
    FirstHex(d,b);
    return c==d;
}

long ProbPrime(const EE& a, long NTRY)
// test if a is eisenstein prime, i.e.,
//   |a|^2 is prime and |a|^2==1 (mod 3) or
//   a.y==0 and |a.x| is prime and |a.x|==2 (mod 3) or
//   a.x==0 and |a.y| is prime and |a.y|==2 (mod 3)
// NTRY = number of trials of Miller-Rabin test
{
    ZZ b;
    if(IsZero(a.y)) abs(b, a.x);
    else if(IsZero(a.x)) abs(b, a.y);
    if(!IsZero(b))
        return b%3 == 2 && ProbPrime(b, NTRY);
    norm(b,a);
    return b==3 || b%3 == 1 && ProbPrime(b, NTRY);
}

void GCD(EE& d, const EE& a, const EE& b)
// d = greatest common divisor of a and b
//   in first hexant 0 < arg(d) < pi/3
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
//   in first hexant 0 < arg(d) < pi/3
//   and compute s,t such that d = s*a + t*b
// by extended euclidean algorithm
{
    long e;
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
    e = FirstHex(d,x);
    rot60e(s,s,e);
    rot60e(t,t,e);
}

void RandomBnd(EE& a, const ZZ& n)
// a = random Eisenstein integer x+yw
//   such that |x|<n and |y|<n
{
    RandomBnd(a.x, n);
    RandomBnd(a.y, n);
    rot60e(a, a, RandomBnd(3)<<1);
}

void GenPrime(EE& p, long l, long f, long err)
// generate random Eisenstein prime
// p = Eisenstein integer such that norm(p)==q^f
//   where q is random prime integer
// l = bit length of q so that 2^{l-1} < q < 2^l
// f = 1 or 2 (q == 1 or 2 mod 3, respectively)
//   If f==1, p is imaginary (and primary),
//   If f==2, p = q (real and positive).
// err = bound for error probability < 2^{-err}.
// If f<1 or f>2, f is set to 1 or 2, respectively.
{
    if(l<2) Error("l<2 in GenPrime");
    ZZ q;
    if(f<=1) {
        do GenPrime(q,l,err);
        while(q%3 == 2);
        if(q==3) { set(p,1,-1); return; }
        FactorPrime(p,q);
        if(RandomBits_long(1)) conj(p,p);
    }
    else {
        do GenPrime(p.x, l, err);
        while(p.x%3 != 2);
        clear(p.y);
    }
}

std::ostream& operator<<(std::ostream& s, const EE& a) {
    s << '[' << a.x << ' ' << a.y << ']';
    return s;
}
// uses NTL
//   http://www.shoup.net/ntl

#include "EEFactoring.h"
using namespace NTL;

void FactorPrime(EE& f, const ZZ& p)
// find x,y such that x^2 - xy + y^2 = p
// where p is prime and p==1 (mod 3)
// return f = x+yw, x==2,y==0(mod 3)
{
    ZZ a,b,c,q,r,t, &x(f.x), &y(f.y);
    sub(b,a=p,3);
    SqrRootMod(b,b,p);
    if(!IsOdd(b)) sub(b,p,b);
    sqr(c,b); c+=3;
    c/=p; c>>=2;
    set(f);
    for(;;) {
        if(a>c) {
            negate(b,b);
            swap(a,c);
            negate(y,y);
            swap(x,y);
        }
        negate(t,a);
        if(b<=a && b>t) break;
        LeftShift(t,a,1);
        DivRem(q,r,b,t);
        if(r>a) { r-=t; q++; }
        add(t,b,r);
        t>>=1;
        b=r;
        MulSubFrom(c,t,q);
        MulAddTo(x,y,q);
    }
    if(sign(b)>0) negate(y,y);
    primary(f,f);
}

void factor(Vec<Pair<EE, long> >& f, const EE& a)
// f = factorization of a into Eisenstein primes
// each element of f is a pair of prime and its exponent
// such that product of prime^{exponent} is associate of a.
// real factors are inserted first in f (if any)
// imaginary factors are appended after real factors
// real factors are positive and sorted in increasing order
// imaginary factors are primary and sorted by norm
{
    int i,j,k(0),e(0);
    ZZ s,t;
    EE b(a);
    Vec<Pair<ZZ, long> > g,h;
    f.SetLength(0);
    if(IsZero(a) || IsUnit(a)) return;
    GCD(t, a.x, a.y);
    b.x /= t;
    b.y /= t;
    norm(s,b);
    factor(g,s);
    factor(h,t);
    for(i=j=0; i<h.length(); i++) {
        if(h[i].a%3 == 2) {
            f.SetLength(k+1);
            f[k].a = h[i].a;
            f[k++].b = h[i].b;
        }
        else { if(i>j) h[j] = h[i]; j++; }
    }
    h.SetLength(j);
    i=j=0;
    if(g.length() && g[0].a == 3) e += g[i++].b;
    if(h.length() && h[0].a == 3) e += h[j++].b*2;
    if(e) {// factor 1-w
        f.SetLength(k+1);
        set(f[k].a, 1, -1);
        f[k++].b = e;
    }
    while(i<g.length() || j<h.length()) {// imaginary factor
        if(j==h.length() || i<g.length() && g[i].a < h[j].a) {
            f.SetLength(k+1);
            FactorPrime(f[k].a, g[i].a);
            if(!divide(b, f[k].a)) conj(f[k].a, f[k].a);
            f[k++].b = g[i++].b;
        }
        else {
            f.SetLength(k+2);
            FactorPrime(f[k].a, h[j].a);
            conj(f[k+1].a, f[k].a);
            f[k+1].b = f[k].b = h[j].b;
            if(i<g.length() && g[i].a == h[j].a) {
                if(divide(b, f[k].a)) f[k].b += g[i++].b;
                else f[k+1].b += g[i++].b;
            }
            j++; k+=2;
        }
    }
}

void mul(EE& a, const Vec<Pair<EE, long> >& f)
// a = product of (Eisenstein integer)^{exponent} in f
// each element of f is a pair of integer and exponent
{
    int i;
    EE b;
    set(a);
    for(i=0; i<f.length(); i++) {
        power(b, f[i].a, f[i].b);
        a *= b;
    }
}
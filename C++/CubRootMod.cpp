// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ_pXFactoring.h>
#include<NTL/ZZ_pEXFactoring.h>
#include "EE.h"

void mod(ZZ_p& b, const EE& a, const EE& p)
// b==a (mod p) where p is eisenstein prime
// Assume norm(p)==1 (mod 3)
// Assume ZZ_p::init(norm(p)) has been executed
{
    ZZ_p x,y;
    conv(x, p.x);
    conv(y, p.y); x /= y;
    conv(y, a.y); y *= x;
    conv(x, a.x); sub(b, x, y);
}

void CubRootMod(EE& x, const EE& a, const EE& p)
// solve x^3 == a (mod p)
// Assume p is primary prime and (a/p)_3 == 1
// Assume either
//   norm(p) is prime and norm(p)==1 (mod 3)
//   or norm(p) is prime^2 and p==2 (mod 3)
{
    ZZ_p c;
    ZZ_pX f;

    if(!IsZero(p.y)) {// norm(p)==1 (mod 3)
        ZZ n;
        norm(n,p);
        ZZ_pPush _p(n);
        mod(c,a,p);
        negate(c,c);
        SetCoeff(f,3);
        SetCoeff(f,0,c);
        FindRoot(c,f);
        conv(n,c);
        conv(x,n);
        x %= p;
    }
    else {// p==2 (mod 3)
        ZZ_pE b;
        ZZ_pEX g;
        ZZ_pPush _p(p.x);
        SetCoeff(f,2);
        SetCoeff(f,1);
        SetCoeff(f,0);
        ZZ_pEPush _f(f);
        clear(f);
        conv(c, a.x); SetCoeff(f,0,c);
        conv(c, a.y); SetCoeff(f,1,c);
        conv(b,f); negate(b,b);
        SetCoeff(g,0,b);
        SetCoeff(g,3);
        FindRoot(b,g);
        conv(f,b);
        conv(x.x, f[0]);
        conv(x.y, f[1]);
    }
}
#include "EEFactoring.h"
#include<iostream>
using namespace NTL;

main() {
    // Greatest Common Divisor
    EE a(2,3),b(4,5),c(6,7);
    EE d,s,t;
    Vec<Pair<EE,long> > f;
    a *= c; b *= c;
    GCD(d,a,b); std::cout << d << std::endl;
    XGCD(d,s,t,a,b);
    s *= a; t *= b; s += t;
    std::cout << (d==s) << std::endl;

    // Factoring integers into Eisenstein primes
    GenPrime(a,8);
    GenPrime(b,8);
    GenPrime(c,8);
    power(a,a,2);
    power(b,b,3); a *= b;
    power(c,c,4); a *= c;
    factor(f,a); std::cout << f << std::endl;
    mul(b,f); std::cout << IsAssoc(a,b) << std::endl;

    // Cubic Residue Symbol
    EE p,x;
    ZZ m;
    GenPrime(p,20);
    GenPrime(a,20); a%=p;
    ResSymb(s,a,p);
    norm(m,p); m--; m/=3;
    PowerMod(t,a,m,p);
    std::cout << s << std::endl;
    std::cout << (s==t) << std::endl;
    if(IsOne(s)) {// solve x^3 == a (mod p)
        CubRootMod(x,a,p);
        std::cout << x << std::endl;
        PowerMod(b,x,3,p);
        std::cout << (a==b) << std::endl;
    }    
}
#ifndef __EE_h__
#define __EE_h__

#include<NTL/ZZ.h>
using namespace NTL;

struct EE {
    // Eisenstein Integer x + yw (w = cubic root of 1)
    ZZ x,y;
    EE() {;}// x=y=0
    EE(const ZZ& a, const ZZ& b) : x(a), y(b) {;}// x=a,y=b
    EE(long a, long b) : x(a), y(b) {;}// x=a,y=b
    EE(const ZZ& a) : x(a) {;}// x=a,y=0
    EE(long a) : x(a) {;}// x=a,y=0
    EE& operator=(const ZZ& a);// x=a
    EE& operator=(long a);// x=a
};

std::ostream& operator<<(std::ostream& s, const EE& a);// output a to s

void conv(EE& b, const ZZ& a);// b = a+0w
void conv(EE& b, long a);// b = a+0w

EE& operator+=(EE& b, const EE& a);// b+=a
EE& operator-=(EE& b, const EE& a);// b-=a
EE& operator*=(EE& b, const EE& a);// b*=a
EE& operator/=(EE& b, const EE& a);// b/=a
EE& operator%=(EE& b, const EE& a);// b%=a

inline long operator==(const EE& a, const EE& b) { return a.x==b.x && a.y==b.y; }
inline long operator!=(const EE& a, const EE& b) { return a.x!=b.x || a.y!=b.y; }
inline long operator==(const EE& a, long b) { return a.x==b && IsZero(a.y); }
inline long operator!=(const EE& a, long b) { return a.x!=b || !IsZero(a.y); }
inline long IsZero(const EE& a) { return IsZero(a.x) && IsZero(a.y); }
inline long IsOne(const EE& a) { return IsOne(a.x) && IsZero(a.y); }

void set(EE& a);// a=1
void set(EE& a, const ZZ& x, const ZZ& y);// a=x+yw
void set(EE& a, long x, long y);// a=x+yw
void clear(EE& a);// a=0
void conj(EE& b, const EE& a);// b = complex conjugate of a
void mirror(EE&, const EE&);// b = reflection of a wrt line y=x/sqrt(3)
void norm(ZZ& x, const EE& a);// b = |a|^2
void negate(EE& b, const EE& a);// b = -a
void rot60(EE& b, const EE& a);// b = a*(1+w)
void rot60(EE& b, const EE& a, long k);// b = a*(1+w)^k
void rot120(EE& b, const EE& a);// b = a*w

long hexant(const EE& a);
// return -1 if a==0
// return k=0,1,...,5 if 60k <= arg(a) < 60(k+1)

long FirstHex(EE& b, const EE& a);
// b = a*(1+w)^k for some k=0,1,...,5
//     in first hexant 0 <= arg(b) < 60
// return k
// if a is zero, then b=0 and k=0

void add(EE& c, const EE& a, const EE& b);// c=a+b
void sub(EE& c, const EE& a, const EE& b);// c=a-b
void mul(EE& c, const EE& a, const EE& b);// c=a*b
void sqr(EE& b, const EE& a);// b=a*a
void div(EE& q, const EE& a, const EE& b);// q=a/b
void rem(EE& r, const EE& a, const EE& b);// r=a%b, norm(r) < norm(b)
void DivRem(EE& q, EE& r, const EE& a, const EE& b);// q=a/b, r=a%b

long divide(EE& q, const EE& a, const EE& b);
// if a/b is divisible, set q=a/b and return 1
// else return 0 and q is unchanged

long divide(const EE& a, const EE& b);
// if a/b is divisible, return 1, else return 0

long divide3(EE& q, const EE& a);
// if a is divisible by 1-w, set q=a/(1-w) and return 1
// else return 0 and q is unchanged

void power(EE& b, const EE& a, long n);// b = a^n (n>=0)

long IsUnit(const EE& a);// test if norm(a)==1

long IsAssoc(const EE& a, const EE& b);
// return 1 if a = b*(1+w)^k for some k=0,1,...,5
// else return 0

void GCD(EE& d, const EE& a, const EE& b);
// d = greatest common divisor of a and b
//     in first hexant 0 < d.y < d.x
// by euclidean algorithm

void XGCD(EE&, EE&, EE&, const EE&, const EE&);
// d = greatest common divisor of a and b
//     in first hexant 0 < d.y < d.x
//     and compute s,t such that d = s*a + t*b
// by extended euclidean algorithm

// b = random Eisenstein integer in hexagon R*e^{i\pi n/3} (n=0...5)
void RandomBits(EE& b, long l);// 0 <= R < 2^l
void RandomLen(EE& b, long l);// 2^{l-1} <= R < 2^l
void RandomBnd(EE& b, const ZZ& a);// 0 <= R < a

void PowerMod(EE& b, const EE& a, long n, const EE& m);
void PowerMod(EE& b, const EE& a, const ZZ& n, const EE& m);
// b = a^n mod m; assume n>=0 and norm(a) < norm(m)

long ProbPrime(const EE& a, long NTRY=10);
// return 1 if either
//   |a| is prime and |a|==2 (mod 3)
//   or norm(a) is prime and norm(a)==1 (mod 3)
//   or norm(a) == 3
// else return 0
// NTRY: number of trials of Miller-Rabin test
// primality is tested probabilistically by MR method

void GenPrime(EE& p, long l, long f=1, long err=80);
// generate random Eisenstein prime p
// Assume f=1 or 2; Assume l>=2
// if f=1, |p|^2=q and q==1 (mod 3)
// if f=2, p=q+0i and q==2 (mod 3)
// where q is random prime and 2^{l-1} <= q < 2^l
// probability of error is less than 2^-err
// p is primary, i.e., p.x==2 and p.y==0 (mod 3)

long primary(EE& b, const EE& a);
// b = unit * a such that b.x==2 and b.y==0 (mod 3)
// return k=0,1,2 such that a = \pm w^k * b
// Assume norm(a) != 0 (mod 3)

void ResSymb(EE& s, const EE& a, const EE& b);
// s = cubic residue symbol (a/b)_3 = 0,1,w,w^2
// Assume norm(a) < norm(b) and norm(b) != 0 (mod 3)
// Assume b is primary, but may not be prime

void CubRootMod(EE&, const EE&, const EE&);
// solve x^3 == a (mod p)
// Assume p is primary prime and (a/p)_3 == 1
// Assume either
//   norm(p) is prime and norm(p)==1 (mod 3)
//   or norm(p) is prime^2 and p==2 (mod 3)

void FactorPrime(EE& f, const ZZ& p);
// find x,y such that x^2 - xy + y^2 = p and x==2, y==0 (mod 3)
// Assume p is prime and p==1 (mod 3)
// return f = x+wy

#endif // __EE_h__
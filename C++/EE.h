// uses NTL
//   http://www.shoup.net/ntl

#ifndef __EE_h__
#define __EE_h__

#include<NTL/ZZ.h>

struct EE {
    // Eisenstein Integer x + yw (w = cubic root of 1)
    NTL::ZZ x,y;
    EE() {;}
    EE(const NTL::ZZ& a, const NTL::ZZ& b) : x(a), y(b) {;}
    EE(long a, long b) : x(a), y(b) {;}
    EE(const NTL::ZZ& a) : x(a) {;}// a+0i
    EE(long a) : x(a) {;}// a+0i
    EE& operator=(const NTL::ZZ& a);// a+0i
    EE& operator=(long a);// a+0i
};

EE& operator+=(EE& b, const EE& a);
EE& operator-=(EE& b, const EE& a);
EE& operator*=(EE& b, const EE& a);
EE& operator/=(EE& b, const EE& a);
EE& operator%=(EE& b, const EE& a);

inline long operator==(const EE& a, const EE& b) { return a.x==b.x && a.y==b.y; }
inline long IsZero(const EE& a) { return NTL::IsZero(a.x) && NTL::IsZero(a.y); }
inline long IsOne(const EE& a) { return NTL::IsOne(a.x) && NTL::IsZero(a.y); }

void set(EE& a);// a=1
void set(EE& a, const NTL::ZZ& x, const NTL::ZZ& y);// a=x+yw
void set(EE& a, long x, long y);// a=x+yw
void clear(EE& a);// a=0
void conj(EE& b, const EE& a);// complex conjugate
void mirror(EE& b, const EE& a);// reflect wrt line y=x/sqrt(3)
void norm(NTL::ZZ& n, const EE& a);// n = a*conj(a)
void negate(EE& b, const EE& a);// b=-a
void rot60(EE& b, const EE& a);// b = a*(1+w)
void rot120(EE& b, const EE& a);// b = a*w
void rot60e(EE& b, const EE& a, long e);// b = a*(1+w)^e

long hexant(const EE& a);
// return floor(arg(a)*3/pi) if a!=0
// return -1 if a==0

long FirstHex(EE& b, const EE& a);
// b = a*(1+w)^e so that 0<=arg(b)<pi/3
// and return e

void primary(EE& b, const EE& a);
// b = a*unit == x+yw to make x==2 and y==0(mod 3)
// assume norm(a) != 0 (mod 3)

void add(EE& c, const EE& a, const EE& b);// c=a+b
void sub(EE& c, const EE& a, const EE& b);// c=a-b
void mul(EE& c, const EE& a, const EE& b);// c=a*b
void sqr(EE& b, const EE& a);// b=a*a

void div(EE& q, const EE& a, const EE& b);
// q = quotient of a/b such that
//   a = bq + r and |r/b|^2 <= 3/4

void rem(EE& r, const EE& a, const EE& b);
// r = remainder of a/b such that
//   a = bq + r and |r/b|^2 <= 3/4

void DivRem(EE& q, EE& r, const EE& a, const EE& b);
// q,r = quotient and remainder of a/b such that
//   a = bq + r and |r/b|^2 <= 3/4

long divide(EE& q, const EE& a, const EE& b);
// if a/b is divisible, set q=a/b and return 1
// else return 0 (and q is unchanged)

long divide(const EE& a, const EE& b);
// if a/b is divisible, return 1, else return 0

void power(EE& b, const EE& a, long e);
// b = a^n; assume n>=0

long IsUnit(const EE& a);// test if norm(a)==1

long IsAssoc(const EE& a, const EE& b);
// return 1 if a = b*w^e for some e (e=0,1,2,3) else 0

long ProbPrime(const EE& a, long NTRY=10);
// test if a is eisenstein prime, i.e.,
//   |a|^2 is prime and |a|^2==1 (mod 3) or
//   a.y==0 and |a.x| is prime and |a.x|==2 (mod 3) or
//   a.x==0 and |a.y| is prime and |a.y|==2 (mod 3)
// NTRY = number of trials of Miller-Rabin test

void GCD(EE&, const EE&, const EE&);
// d = greatest common divisor of a and b
//   in first hexant 0 < arg(d) < pi/3
// by euclidean algorithm

void XGCD(EE&, EE&, EE&, const EE&, const EE&);
// d = greatest common divisor of a and b
//   in first hexant 0 < arg(d) < pi/3
//   and compute s,t such that d = s*a + t*b
// by extended euclidean algorithm

void RandomBnd(EE& a, const NTL::ZZ& n);
// a = random Eisenstein integer x+yw
//   such that |x|<n and |y|<n

void GenPrime(EE& p, long l, long f=1, long=80);
// generate random Eisenstein prime
// p = Eisenstein integer such that norm(p)==q^f
//   where q is random prime integer
// l = bit length of q so that 2^{l-1} < q < 2^l
// f = 1 or 2 (q == 1 or 2 mod 3, respectively)
//   If f==1, p is imaginary (and primary),
//   If f==2, p = q (real and positive).
// err = bound for error probability < 2^{-err}.
// If f<1 or f>2, f is set to 1 or 2, respectively.

std::ostream& operator<<(std::ostream& s, const EE& a);
// print a as [a.x a.y]

void FactorPrime(EE&, const NTL::ZZ&);
// find x,y such that x^2 - xy + y^2 = p
// where p is prime and p==1 (mod 3)
// return f = x+yw, y==0 (mod 3)

#endif // __EE_h__
// uses NTL
//   http://www.shoup.net/ntl

#ifndef __EEFactoring_h__
#define __EEFactoring_h__

#include<NTL/pair.h>
#include "EE.h"

void factor(NTL::Vec<NTL::Pair<NTL::ZZ, long> >& f, const NTL::ZZ& n);
// find x,y such that x^2 - xy + y^2 = p
// where p is prime and p==1 (mod 3)
// return f = x+yw, x==2,y==0(mod 3)

void factor(NTL::Vec<NTL::Pair<EE, long> >& f, const EE& a);
// f = factorization of a into Eisenstein primes
// each element of f is a pair of prime and its exponent
// such that product of prime^{exponent} is associate of a.
// real factors are inserted first in f (if any)
// imaginary factors are appended after real factors
// real factors are positive and sorted in increasing order
// imaginary factors are primary and sorted by norm

void mul(NTL::ZZ& a, const NTL::Vec<NTL::Pair<NTL::ZZ, long> >& f);
// a = product of (integer)^{exponent} in f
// each element of f is a pair of integer and exponent

void mul(EE& a, const NTL::Vec<NTL::Pair<EE, long> >& f);
// a = product of (Eisenstein integer)^{exponent} in f
// each element of f is a pair of integer and exponent

#endif // __GGFactoring_h__

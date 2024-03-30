from EE import *

a,b,c = EE(2,3),EE(4,5),EE(6,7)
a *= c; b *= c
d = GCD(a,b); print(d)
d,s,t = XGCD(a,b)
print(d == s*a + t*b)

a,b,c = GenPrime(8),GenPrime(8),GenPrime(8)
a = a**2 * b**3 * c**4; print(a)
f = factor(a); print(f)
b = 1
for k in f: b *= k**f[k]
print(IsAssoc(a,b))

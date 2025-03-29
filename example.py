from EE import *

# Greatest Common Divisor
a,b,c = EE(2,3),EE(4,5),EE(6,7)
a *= c; b *= c
d = GCD(a,b); print(d)
d,s,t = XGCD(a,b)
print(d == s*a + t*b)

# Factoring integers into Eisenstein primes
a,b,c = GenPrime(8),GenPrime(8),GenPrime(8)
a = a**2 * b**3 * c**4; print(a)
f = factor(a); print(f)
b = 1
for k in f: b *= k**f[k]
print(IsAssoc(a,b))

# Cubic Residue Symbol
p = GenPrime(20)
a = GenPrime(20); a %= p
s = ResSymb(a,p); print(a,p,s)
print(s == PowerMod(a, (norm(p)-1)//3, p))
if s==1:# solve x^3 == a (mod p)
    x = CubRootMod(a,p); print(x)
    print(PowerMod(x,3,p) == a)


from math import sqrt
w = (-1 + sqrt(3)*1j)/2 # cube root of unity

class EE:
    """ Eisenstein integer x+yw """
    def __init__(a, x=0, y=0):
        a.x = int(x)
        a.y = int(y)

    def __repr__(a):
        if a.y == 0: return str(a.x)
        elif a.x == 0: return str(a.y) + 'w'
        elif a.y < 0:
            return '(' + str(a.x) + str(a.y) + 'w)'
        else:
            return '(' + str(a.x) + '+' + str(a.y) + 'w)'

    def __hash__(a):
        return a.x ^ a.y

    def __neg__(a):
        return EE(-a.x, -a.y)

    def __add__(a,b):
        if isinstance(b, EE):
            return EE(a.x + b.x, a.y + b.y)
        else:
            return EE(a.x + b, a.y)

    def __sub__(a,b):
        if isinstance(b, EE):
            return EE(a.x - b.x, a.y - b.y)
        else:
            return EE(a.x - b, a.y)

    def __mul__(a,b):
        if isinstance(b, EE):
            s = a.x*b.x
            t = a.y*b.y
            u = (a.x - a.y)*(b.y - b.x)
            return EE(s-t, s+u)
        else:
            return EE(b*a.x, b*a.y)

    def __truediv__(a,b):
        """ a = bq + r, norm(r) <= norm(b)/2 """
        if isinstance(b, EE):
            a *= conj(b)
            b = norm(b)
        x = ((a.x<<1) + b)//(b<<1)
        y = ((a.y<<1) + b)//(b<<1)
        return EE(x,y)

    def __mod__(a,b):
        q = a/b
        return a - b*q

    def __divmod__(a,b):
        q = a/b
        return q, a - b*q

    def __pow__(a,k): # assume k>=0
        if k==0: return EE(1)
        n = (1<<(k.bit_length()-1))>>1
        b = a
        while n:
            b *= b
            if k&n: b *= a
            n>>=1
        return b

    def __eq__(a,b):
        if isinstance(b, EE):
            return a.x == b.x and a.y == b.y
        else:
            return a.x == b and a.y == 0

    def __ne__(a,b):
        return not a.__eq__(b)

    def __radd__(a,b):
        return EE(b + a.x, a.y)

    def __rsub__(a,b):
        return EE(b - a.x, -a.y)

    def __rmul__(a,b):
        return EE(b*a.x, b*a.y)

    def __rtruediv__(a,b):
        return EE(b).__truediv__(a)

    def __rmod__(a,b):
        return EE(b).__mod__(a)

    def __rdivmod__(a,b):
        return EE(b).__divmod__(a)

    def __complex__(a):
        return a.x + w*a.y

    def isUnit(a):# test if norm(a)==1
        return (a.y==0 and abs(a.x)==1) or\
               (abs(a.y)==1 and (a.x==0 or a.x==a.y))

def norm(a): return (a.x - a.y)**2 + a.x*a.y# |a|^2
def conj(a): return EE(a.x-a.y, -a.y) # complex conjugate of a
def rot60(a): return EE(a.x-a.y, a.x) # a*(1+w)
def rot120(a): return EE(-a.y, a.x-a.y) # a*w
def mirror(a): return EE(a.x, a.x-a.y) # reflect wrt y=x/sqrt(3)

def IsAssoc(a, b, exponent=False):
    """ a,b: EE, return int
    test if a == b*(1+w)^k for some k=0,1,...,5
    if exponent is False, return True or False
    else if a is associate of b, return k
    else return -1
    """
    c = FirstHex(a, exponent)
    d = FirstHex(b, exponent)
    if exponent:
        if c[0]!=d[0]: return -1
        k = d[1] - c[1]
        return k if k>=0 else k+6
    else: return c==d

def rot60k(a,k):# a*(1+w)^k
    """ a: EE, k: int, return EE """
    k %= 6
    if   k==1: return rot60(a)
    elif k==2: return rot120(a)
    elif k==3: return -a
    elif k==4: return -rot60(a)
    elif k==5: return -rot120(a)
    else: return a

def hexant(a):
    """ a: EE, return int
    return -1 if a==0,
    return k if 60k <= arg(a) < 60*(k+1) (k=0,1,...,5)
    """
    if a.y==0:
        if a.x==0: return -1
        elif a.x>0: return 0
        else: return 3
    elif a.y>0:
        if a.x > a.y: return 0
        elif a.x > 0: return 1
        else: return 2
    elif a.x < a.y: return 3
    elif a.x < 0: return 4
    else: return 5

def FirstHex(a, exponent=False):
    """ a: EE, return EE
    return b = a*(1+w)^k for some k=0,1,...,5
      in first hexant 0 <= arg(b) < 60
    if exponent is True, return b and k
    if a is zero, then b=0 and k=0
    """
    k = hexant(a)
    k = 6-k if k>0 else 0
    b = rot60k(a,k)
    if exponent: return b,k
    else: return b

def GCD(a,b):
    """ a,b: EE, return EE
    d = greatest common divisor of a and b
      in first hexant 0 <= arg(d) < 60
      by Euclidean algorithm
    return d
    if a and b are both zero, then d=0
    """
    while b!=0: a,b = b,a%b
    return FirstHex(a)

def XGCD(a,b):
    """ a,b: EE, return EE*3
    d = greatest common divisor of a and b
      in first hexant 0 <= arg(d) < 60
      by extended Euclidean algorithm
    return d,s,t such that d = s*a + t*b
    if a and b are both zero, then d,s,t=0,1,0
    """
    s,t = EE(1),EE(0)
    u,v = EE(0),EE(1)
    while b!=0:
        q,r = divmod(a,b)
        s,u = u, s-u*q
        t,v = v, t-v*q
        a,b = b,r
    d,k = FirstHex(a, True)
    s = rot60k(s,k)
    t = rot60k(t,k)
    return d,s,t

def primary(a, exponent=False):
    """ a: EE, return EE
    return b such that a = \pm w^k * b (k=0,1,2)
    and b.x==2 (mod 3) and b.y==0 (mod 3)
    Assume norm(a) != 0 (mod 3)
    if exponent is True, return b and k
    """
    x,y = a.x%3, a.y%3
    if y==2: a = -a
    if x==y: a,e = rot120(a), 2
    elif x==0: a,e = rot60(a), 1
    elif x==1: a,e = -a,0
    else: e = 0
    if exponent: return a,e
    else: return a

def ResSymb(a,b):
    """a,b: EE, return EE
    return cubic residue symbol (a/b)_3 = 0,1,w,w^2
    Assume norm(a) < norm(b) and norm(b) != 0 (mod 3)
    Assume b is primary, but may not be prime
    reference: K. Ireland and M. Rosen
      "A Classical Introduction to Modern Number Theory" section 9.3
    """
    j,lam = 0, EE(1,-1)
    while a!=0:
        m,n = (b.x+1)//3, b.y//3
        m,n = m%3,(m+n)%3
        while True:
            q,r = divmod(a, lam)
            if r!=0: break
            a,j = q, j-m

        a,k = primary(a, True)
        a,b,j = b%a, a, (j+k*n)%3

    if not b.isUnit(): return EE()
    elif j==0: return EE(1)
    elif j==1: return EE(0,1)
    else: return EE(-1,-1)

def PowerMod(a,k,n):# a^k mod n
    """ a,n: EE, k: int, return EE """
    if k==0: return EE(1)
    m = (1<<(k.bit_length()-1))>>1
    b = a
    while m:
        b = b*b%n
        if k&m: b = b*a%n
        m>>=1
    return b

#########################################################
import sympy as sp
from sympy.abc import x
from math import gcd
from random import randrange

def FactorPrime(p):
    """ p: int, return EE
    given prime number p,
    find x,y such that x*x - x*y + y*y = p
         and x==2, y==0 (mod 3)
    Assume p==1 (mod 3)
    return f = x+y*w (primary prime)
    """
    a,b = p, sp.sqrt_mod(p-3, p)
    if b&1==0: b = p-b
    c,x,y = ((b*b + 3)//p)>>2, 1,0
    while True:
        if a>c: a,b,c,x,y = c,-b,a,-y,x
        if b<=a and b>-a: break
        q,r = divmod(a-b, a<<1)
        b,c,x = a-r, c+((a+b-r)>>1)*q, x-q*y

    return primary(EE(x,-y))

def factor(a):
    """ a: EE, return dict{EE,int}
    factorize a into Eisenstein primes
    return dictionary of (prime, exponent) pair
    such that product of p**e is associate of a
    real factors are positive and
    imaginary factors are in primary
    """
    if not isinstance(a,EE): a = EE(a)
    f = {}
    if a==0 or a.isUnit(): return f
    d = gcd(a.x, a.y)
    if d>1: a /= d
    g = sp.factorint(norm(a))
    h = sp.factorint(d)

    k = g.pop(3,0) + 2*h.pop(3,0)
    if k: f[EE(1,-1)] = k
    for k in h:
        if k%3 == 2: f[k] = h[k]
        else:
            p = FactorPrime(k)
            q = conj(p)
            f[p],f[q] = h[k],h[k]
            if k in g:
                if a%p == 0: f[p] += g.pop(k)
                else:        f[q] += g.pop(k)
    for k in g:
        p = FactorPrime(k)
        if a%p != 0: p = conj(p)
        f[p] = g[k]

    return f

def IsPrime(a):
    """ a: EE, return bool
    test if a is Eisenstein prime
    """
    if not isinstance(a,EE): a = EE(a)
    if a.y==0 or a.x==a.y: b = abs(a.x)
    elif a.x==0: b = abs(a.y)
    else: b = -1
    if b>=0: return b%3 == 2 and sp.isprime(b)
    b = norm(a)
    return b==3 or (b%3 == 1 and sp.isprime(b))

def GenPrime(l, f=1):
    """ l,f: int, return EE
    generate random Eisenstein prime p
    Assume f=1 or 2; Assume l>=2
    if f=1, |p|^2=q and q==1 (mod 3)
    if f=2, p=q+0i and q==2 (mod 3)
    where q is random prime and 2^{l-1} <= q < 2^l
    p is primary, i.e., p.x==2, p.y==0 (mod 3)
    """
    while True:
        p = sp.randprime(1<<(l-1), 1<<l)
        if f <= 1:
            if p == 3: return EE(1,-1)
            if p%3 == 1:
                p = FactorPrime(p)
                if randrange(2): return p
                else: return conj(p)

        elif p%3 == 2: return EE(p)

def InvMod(a,p):
    """ a,p: int, return int
    return x such that ax==1 (mod p)
    Assume a and p are relatively prime
    """
    s,u = 1,0
    while p:
        q,r = divmod(a,p)
        a,p,s,u = p,r,u,s-q*u
    return s

def CubRootMod(a,p):
    """ a,p: EE, return EE
    solve x^3 == a (mod p) and return x
    Assume norm(p) is prime and norm(p)==1 (mod 3)
    Assume (a/p)_3 == 1
    """
    n = norm(p)
    b = p.x * InvMod(p.y, n) % n
    b = (a.x - b * a.y % n) % n
    F = sp.factor_list(x**3 - b, modulus=n)
    return int(-F[1][0][0].coeff(x,0))%p

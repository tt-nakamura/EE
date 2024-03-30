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
        """ a = bq + r, norm(r) <= (3/4)*norm(b) """
        if isinstance(b, EE):
            a *= conj(b)
            b = norm(b)

        x = ((a.x<<1) + b)//(b<<1)
        y = ((a.y<<1) + b)//(b<<1)
        return EE(x,y)

    def __mod__(a,b):
        return a - a/b*b

    def __divmod__(a,b):
        q = a/b
        return q, a - b*q

    def __pow__(a,e):
        """ assume e>=0 """
        if e==0: return EE(1)
        n = (1<<(e.bit_length()-1))>>1
        b = a
        while n:
            b *= b
            if e&n: b *= a
            n>>=1
        return b

    def __eq__(a,b):
        if isinstance(b, EE):
            return a.x == b.x and a.y == b.y
        else:
            return a.x == b and a.y == 0

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

def norm(a): return (a.x - a.y)**2 + a.x*a.y

def conj(a): return EE(a.x-a.y, -a.y)
def IsUnit(a): return (a.y==0 and abs(a.x)==1)\
                   or (a.y==1 and (a.x==0 or a.x==1))\
                   or (a.y==-1 and (a.x==0 or a.x==-1))
def rot60(a): return EE(a.x-a.y, a.x) # -a/w
def rot120(a): return EE(-a.y, a.x-a.y) # a*w
def mirror(a): return EE(a.x, a.x-a.y) # reflect wrt y=x/sqrt(3)

def IsAssoc(a, b, exponent=False):
    """ test if a == b*w^{e/2} for some e (e=0,1,...,5)
    if exponent is False, return True or False
    else if a is associate of b, return e
    else return -1
    """
    c = FirstHex(a, exponent)
    d = FirstHex(b, exponent)
    if exponent:
        if c[0]!=d[0]: return -1
        e = d[1] - c[1]
        return e if e>=0 else e+6
    else: return c==d

def rot60e(a,e):
    """ a*w^(e/2) """
    e %= 6
    if   e==1: return rot60(a)
    elif e==2: return rot120(a)
    elif e==3: return -a
    elif e==4: return -rot60(a)
    elif e==5: return -rot120(a)
    else: return a

def hexant(a):
    """ return -1 if a==0,
    return e if 60e <= arg(a) < 60*(e+1) (e=0,1,...,5)
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
    """ return b = a*w^(e/2) for some e (e=0,1,...,5)
    in first hexant 0 <= arg(b) < 60
    if exponent is True, return b and e
    if a is zero, then b=0 and e=0
    """
    e = hexant(a)
    e = 6-e if e>0 else 0
    b = rot60e(a,e)
    if exponent: return b,e
    else: return b

def GCD(a,b):
    """ d = greatest common divisor of a and b
    in first hexant 0 <= arg(d) < 60
    by Euclidean algorithm
    if a and b are both zero, then d=0
    """
    while b!=0: a,b = b,a%b
    return FirstHex(a)

def XGCD(a,b):
    """ d = greatest common divisor of a and b
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
    d,e = FirstHex(a, True)
    s = rot60e(s,e)
    t = rot60e(t,e)
    return d,s,t

#########################################################
from sympy import sqrt_mod, factorint, randprime, isprime
from math import gcd
from random import randrange

def primary(a):
    """ return a*unit = x+yw
      such that x==2 and y==0 (mod 3).
    assume norm(a) != 0 (mod 3)
    """
    x,y = a.x%3, a.y%3
    if y==2: a = -a
    if x==y: return rot120(a)
    elif x==0: return rot60(a)
    elif x==1: return -a
    else: return a

def FactorPrime(p):
    """ given a prime number p where p==1 (mod 3),
    find x,y such that x*x - x*y + y*y = p
    and return primary(x+yw)
    """
    a,b = p, sqrt_mod(p-3, p)
    if b&1==0: b = p-b
    c,x,y = ((b**2 + 3)//a)>>2, 1,0
    while True:
        if a>c: a,b,c,x,y = c,-b,a,-y,x
        if b<=a and b>-a: break
        q,r = divmod(a-b, a<<1)
        b,c,x = a-r, c+((a+b-r)>>1)*q, x-q*y

    return primary(EE(x,-y))

def factor(a):
    """ factorize a into Eisenstein primes
    return dictionary of (prime, exponent) pair
    such that product of p**e is associate of a
    real factors are positive and
    imaginary factors are primary
    """
    if not isinstance(a,EE): a = EE(a)
    f = {}
    if a==0 or IsUnit(a): return f
    d = gcd(a.x, a.y)
    if d>1: a /= d
    g = factorint(norm(a))
    h = factorint(d)

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
    """ test if a is Eisenstein prime """
    if not isinstance(a,EE): a = EE(a)
    if   a.y==0: b = abs(a.x)
    elif a.x==0: b = abs(a.y)
    else: b = -1
    if b>=0: return b%3 == 2 and isprime(b)
    b = norm(a)
    return b==3 or (b%3 == 1 and isprime(b))

def GenPrime(l, f=1):
    """ generate random Eisenstein prime p
    such that norm(p) == q^f (f=1,2) where
    q is prime number, bit length of norm(p) is l
    and p is primary; assume l>=2
    """
    while True:
        p = randprime(1<<(l-1), 1<<l)
        if f<=1:
            if p == 3: return EE(1,-1)
            if p%3 == 1:
                p = FactorPrime(p)
                if randrange(2): return p
                else: return conj(p)
        elif p%3 == 2: return p

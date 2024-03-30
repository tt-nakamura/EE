// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_GF2.h>
using namespace NTL;

#define MPQS_MAXLEN 180
#define MPQS_BOUND  0.5
#define MPQS_INTVL  1
#define MPQS_SIEV   1
#define MPQS_EXTRA  10

long Jacobi(long, long);
long SqrRootMod(long, long);

long mpqs(ZZ& d, const ZZ& n)
// input:
//   n = odd integer, not prime power, n>2000
// output:
//   d = divisor of n, 1 < d < n
//       by quadratic sieve method
// return:
//   0 if successful, -1 or -2 if failure
// reference:
//   R. Crandall and C. Pomerance
//     "Prime Numbers: A Computational Perspective"
//     2nd edition (Springer) section 6.1
{
    if(NumBits(n) > MPQS_MAXLEN) return -1;
    long i,j,k,l,m,p,r,s,t,K,Q,B,M,U,N,T;
    double lnN, lnB;
    static double LN2R(1./log(2));
    ZZ a,b,c,q,u;
    ZZ_p x,y,z;
    Vec<long> F,S,sv,rl,FA;
    Vec<char> LF;
    vec_ZZ_p FZ,ru;
    Mat<long> e;
    mat_GF2 A,X;

    if(&d==&n) return mpqs(d,a=n);
    lnN = log(n);
    lnB = MPQS_BOUND*sqrt(lnN*log(lnN));
    B = long(exp(lnB));

    PrimeSeq ps;
    F.append(ps.next());
    S.append(1);
    while((p=ps.next()) <= B) {
        l = n%p;
        if((j = Jacobi(l,p)) > 0) {
            F.append(p);
            S.append(SqrRootMod(l,p));
        }
        else if(j==0) { d=p; return 0; }
    }
    K = F.length();
    M = long(MPQS_INTVL*B);
    U = (M<<1)+1;
    N = K + MPQS_EXTRA;
    ZZ_p::init(n);
    sv.SetLength(U);
    LF.SetLength(K);
    FZ.SetLength(K);
    ru.SetLength(N);
    rl.SetLength(N);
    e.SetDims(N,K+1);

    for(i=0; i<K; i++) LF[i] = char(round(log(F[i])*LN2R));
    for(i=0; i<K; i++) conv(FZ[i], F[i]);
    T = long((0.5*lnN + lnB)*LN2R - MPQS_SIEV*LF[K-1]);
    LeftShift(q,n,1);
    SqrRoot(q,q); q/=M;
    SqrRoot(q,q);
    for(k=0, l=K; k<N; q++) {
        NextPrime(q,q);
        if((j = Jacobi(n,q)) < 0) continue;
        else if(j==0) { d=q; return 0; }
        sqr(a,q); rem(d,n,q);
        SqrRootMod(b,d,q);
        AddMod(c,b,b,q);
        InvMod(d,c,q);
        sqr(c,b); sub(c,n,c); c/=q;
        MulMod(c,c,d,q);
        MulAddTo(b,c,q);
        RightShift(d,a,1);
        if(b>d) sub(b,a,b);
        sqr(c,b); c-=n; c/=a;
        for(i=0; i<U; i++) sv[i] = 0;
        for(j=0; j<K; j++) {
            if(q==(p=F[j])) continue;
            r = InvMod(a%p,p);
            m = b%p;
            for(t=0, s=S[j]; t<2 && p>=3; t++, s=p-s) {
                if((i=((s-m)*r+M)%p)<0) i+=p;
                for(; i<U; i+=p) sv[i] += LF[j];
            }
        }
        for(i=0, s=-M, r=k; i<U && k<N; i++, s++) {
            if(sv[i] < T) continue;
            mul(u,a,s); u+=b;
            add(d,u,b); d*=s; d+=c;
            if(IsZero(d)) continue;
            for(j=0; j<=K; j++) e[k][j] = 0;
            if(sign(d) < 0) e[k][K] = 1;
            abs(d,d);
            for(j=0; j<K; j++) {
                if((p = F[j]) > d) break;
                while(divide(d,d,p)) e[k][j]++;
            }
            if(!IsOne(d)) continue;
            conv(ru[k], u);
            rl[k] = l;
            k++;
        }
        if(k>r) {
            FZ.SetLength(l+1);
            conv(FZ[l++], q);
        }
    }
    F.kill();
    S.kill();
    A.SetDims(N,K+1);
    FA.SetLength(l);
    for(i=0; i<N; i++)
        for(j=0; j<=K; j++) if(e[i][j] & 1) set(A[i][j]);
    kernel(X,A);
    for(k=0; k<X.NumRows(); k++) {
        for(i=0; i<l; i++) FA[i] = 0;
        set(x);
        set(y);
        for(i=0; i<N; i++) {
            if(IsZero(X[k][i])) continue;
            x *= ru[i];
            FA[rl[i]] += 2;
            for(j=0; j<K; j++) FA[j] += e[i][j];
        }
        for(i=0; i<l; i++) {
            if(FA[i] == 0) continue;
            power(z, FZ[i], FA[i]>>1);
            y *= z;
        }
        conv(a,x);
        conv(b,y);
        a -= b;
        GCD(d,a,n);
        if(d>1 && d<n) return 0;
    }
    return -2;
}

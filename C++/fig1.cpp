#include "EE.h"
#include<fstream>

void ResSymb_(EE& s, const EE& a, const EE& p) {
    ZZ n;
    norm(n,p); n--; n/=3;
    PowerMod(s,a,n,p);
}

main() {
    long i,j,k,l,n(100),M(10),N(10),MN(M*N);
    double l1(100),l2(1000);
    std::ofstream f("fig1.txt");
    double t1,t2,s,dl(pow(l2/l1, 1./n));
    EE p,a,b;
    for(i=0; i<=n; i++) {
        l = (long)round(l1 * pow(dl,i));
        t1 = t2 = 0;
        for(j=0; j<M; j++) {
            GenPrime(p,l);
            for(k=0; k<N; k++) {
                RandomLen(a,l); a %= p;
                s = GetTime(); ResSymb_(b,a,p);
                t1 += GetTime() - s;
                s = GetTime(); ResSymb(a,a,p);
                t2 += GetTime() - s;
                if(a!=b) std::cout << "a!=b" << std::endl;
            }
        }
        t1 /= MN;
        t2 /= MN;
        f << l << ' ' << t1 << ' ' << t2 << std::endl;
        std::cout << l << ' ' << t1 << ' ' << t2 << std::endl;
    }
}

/* 
 * File:   randnumbers.h
 * Author: peischl
 *
 * Created on 6. Februar 2013, 18:44
 */

#ifndef RANDNUMBERS_H
#define	RANDNUMBERS_H

/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */
//------------------------------------------------------------------------------
#define MBIG 1000000000L
//Loro_10_10_12 correction suggested by http://www.shadlen.org/ichbin/random/generators.htm
//#define MBIG 2147483647
#define MSEED 161803398L
#define MZ 0
#define FAC (1.0/MBIG)


using namespace std;

double ran3(long *idum)
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;

    if ((*idum < 0) || (iff == 0) ) {
        iff=1;
        mj=MSEED-(*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;++i) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;++k)
            for (i=1;i<=55;++i) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext=0;
        inextp=31;
        *idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;

    //for debug
    if( (mj*FAC) >= 1.0){
      mj=(MBIG-1);
      return mj*FAC;
    }
   else{
      return mj*FAC;
   }
}

//------------------------------------------------------------------------------

void initializeRan3(const long &lidum) {
   long curSeed=-lidum;
   ran3(&curSeed);
   return;
}

//------------------------------------------------------------------------------

#define PI 3.141592654

//double poidev(double xm, long *idum)
//{
//    double gammln(double xx);
//
//    static double sq,alxm,g,oldm=(-1.0);
//    double em,t,y;
//
//    if (xm < 12.0) {
//        if (xm != oldm) {
//            oldm=xm;
//            g=exp(-xm);
//        }
//        em = -1;
//        t=1.0;
//        do {
//            ++em;
//            t *= ran3(idum);
//        } while (t > g);
//    } else {
//        if (xm != oldm) {
//            oldm=xm;
//            sq=sqrt(2.0*xm);
//            alxm=log(xm);
//            g=xm*alxm-gammln(xm+1.0);
//        }
//        do {
//            do {
//                y=tan(PI*ran3(idum));
//                em=sq*y+xm;
//            } while (em < 0.0);
//            em=floor(em);
//            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
//        } while (ran3(idum) > t);
//    }
//    return em;
//}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */



double ran3(long *idum);

void initializeRan3(const long &lidum);

#endif	/* RANDNUMBERS_H */


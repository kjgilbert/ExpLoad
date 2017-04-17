#include "rng2.h"
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;

//------------------------------------------------------------------------------
/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */
#define MBIG 1000000000L
//Loro_10_10_12 correction suggested by http://www.shadlen.org/ichbin/random/generators.htm
//#define MBIG 2147483647
#define MSEED 161803398L
#define MZ 0
#define FAC (1.0/MBIG)

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
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

//------------------------------------------------------------------------------

void initializeRan3(const long &lidum) {
   long curSeed=-lidum;
   ran3(&curSeed);
   return;
}
/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */
//------------------------------------------------------------------------------
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
//------------------------------------------------------------------------------
#define PI 3.141592654
double poidev(double xm, long *idum)
{
	double gammln(double xx);

	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran3(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran3(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran3(idum) > t);
	}
	return em;
}
#undef PI
//------------------------------------------------------------------------------
// randint(a, b) returns a uniformly distributed random integer on [a, b]
int randint(const int & lower_bound, const int & upper_bound) {
   long lidum=1L;
   return int(lower_bound+ran3(&lidum)*((upper_bound-lower_bound)+1));
}
//------------------------------------------------------------------------------
// Returns a uniformly distributed random deviate.  
//	randreal() is uniform on [0, 1[
double randreal() {
   long lidum=1L;
   return ran3(&lidum);
}
//------------------------------------------------------------------------------
//	randreal(a, b) is uniform on [a, b[.
double randreal(const double & lower_bound, const double & upper_bound) {
   long lidum=1L;
   return lower_bound+ran3(&lidum)*(upper_bound-lower_bound);
}
//------------------------------------------------------------------------------
// Returns a random deviate from an exponential distribution with parameter lambda
double randexp(const double & lambda) {
   long lidum=1L;
   return -log(ran3(&lidum))/lambda; 
}
//------------------------------------------------------------------------------
// Returns a random deviate from an exponential distribution with parameter lambda
double randpois(const double & lambda) {
   long lidum=1L;
   return poidev(lambda, &lidum);
}
//------------------------------------------------------------------------------
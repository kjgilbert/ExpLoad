#ifndef __RNG2__
#define __RNG2__

#include <math.h>

using namespace std;

typedef double my_float;
// A routine for generating good random numbers
extern double ran3(long *idum);
void initializeRan3(const long &lidum);

extern double gammln(double xx);
// A routine to generate poisson deviates with mean xm
extern double poidev(double xm, long *idum);

// randint(a, b) returns a uniformly distributed random integer on [a, b]
extern int randint(const int & lower_bound, const int & upper_bound);

// Returns a uniformly distributed random deviate.  
//	randreal() is uniform on [0, 1[
extern double randreal();
//	randreal(a, b) is uniform on [a, b[.
extern double randreal(const double & lower_bound, const double & upper_bound);

// Returns a random deviate from an exponential distribution with parameter lambda
extern double randexp(const double & lambda);

// Returns a random deviate from an exponential distribution with parameter lambda
extern double randpois(const double & lambda);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#endif //__RNG2__
//! \file Gamma function extracted from ROOT
/*! Added for computation of shower profile in CalClustersAlg
 *  \author Regis Terrier
 *  \b Revisions:
 *  -  11/01/00  RT  since abs(double) may not be used in VC++ switch to fabs()
 *  -  05/00    RT    first implementation
 */

#include "gamma.h"


#include <cmath>

  /*! Computation of ln[gamma(z)] for all z>0.
   
    The algorithm is based on the article by C.Lanczos [1] as denoted in
    Numerical Recipes 2nd ed. on p. 207 (W.H.Press et al.).
   
    [1] C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
   
    The accuracy of the result is better than 2e-10.
   
   --- Nve 14-nov-1998 UU-SAP Utrecht
*/
double LnGamma(double z)
{

   if (z<=0) return 0;

   // Coefficients for the series expansion

   double c[7]= { 2.5066282746310005, 76.18009172947146, -86.50532032941677
                   ,24.01409824083091,  -1.231739572450155, 0.1208650973866179e-2
                   ,-0.5395239384953e-5};

   double x   = z;
   double y   = x;
   double tmp = x+5.5;
   tmp = (x+0.5)*log(tmp)-tmp;
   double ser = 1.000000000190015;
   for (int i=1; i<7; i++) {
      y   += 1;
      ser += c[i]/y;
   }
   double v = tmp+log(c[0]*ser/x);
   return v;
}


//_____________________________________________________________________________

  /*! Computation of the incomplete gamma function P(a,x)
  
   The algorithm is based on the formulas and code as denoted in
   Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
  
   --- Nve 14-nov-1998 UU-SAP Utrecht
 */

 double Gamma(double a,double x)
{

   if (a <= 0 || x <= 0) return 0;

   if (x < (a+1)) return GamSer(a,x);
   else           return GamCf(a,x);
}

//______________________________________________________________________________

 /*! Computation of the incomplete gamma function P(a,x)
    via its continued fraction representation.
   
    The algorithm is based on the formulas and code as denoted in
    Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
   
   --- Nve 14-nov-1998 UU-SAP Utrecht
	*/

 double GamCf(double a,double x)
{

   int itmax    = 100;      // Maximum number of iterations
   double eps   = 3.e-7;    // Relative accuracy
   double fpmin = 1.e-30;   // Smallest double value allowed here

   if (a <= 0 || x <= 0) return 0;

   double gln = LnGamma(a);
   double b   = x+1-a;
   double c   = 1/fpmin;
   double d   = 1/b;
   double h   = d;
   double an,del;
   for (int i=1; i<=itmax; i++) {
      an = double(-i)*(double(i)-a);
      b += 2;
      d  = an*d+b;
      if (fabs(d) < fpmin) d = fpmin;
      c = b+an/c;
      if (fabs(c) < fpmin) c = fpmin;
      d   = 1/d;
      del = d*c;
      h   = h*del;
      if (fabs(del-1) < eps) break;
      //if (i==itmax) cout << "*GamCf(a,x)* a too large or itmax too small" << endl;
   }
   
   double v = exp(-x+a*log(x)-gln)*h;
   return (1-v);
}

//______________________________________________________________________________

  /*! Computation of the incomplete gamma function P(a,x)
    via its series representation.
   
    The algorithm is based on the formulas and code as denoted in
    Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
   
   --- Nve 14-nov-1998 UU-SAP Utrecht
 */

 double GamSer(double a,double x)
{
   int itmax  = 100;   // Maximum number of iterations
   double eps = 3.e-7; // Relative accuracy

   if (a <= 0 || x <= 0) return 0;

   double gln = LnGamma(a);
   double ap  = a;
   double sum = 1/a;
   double del = sum;
   for (int n=1; n<=itmax; n++) {
      ap  += 1;
      del  = del*x/ap;
      sum += del;
      if (fabs(del) < fabs(sum*eps)) break;
      //if (n==itmax) cout << "*GamSer(a,x)* a too large or itmax too small" << endl;
   }
   double v = sum*exp(-x+a*log(x)-gln);
   return v;
}

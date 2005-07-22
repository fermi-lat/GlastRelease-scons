/**@file ClassifyCore.cxx
@brief 

*/

#include "ClassifyCore.h"

namespace {
    inline double sqr(double x){return x*x;}
}
double ClassifyCore::scaleFactor(double energy, bool thin)
{
    // following numbers determined empirically to roughly 
    // give a fit scale factor of 1.0 independent of energy
    double t = pow( energy/100., -0.81);
	if( thin )  return sqrt( sqr(28e-3*t) + sqr( 158E-6) ); 
    else        return sqrt( sqr(46e-3*t) + sqr( 358e-6) ); 
}

// $Id$
//

#include "geometry/RootFinder.h"
#include <cmath>
#include <cfloat>

// the static constants
float RootFinder::xacc = (float)1e-5;
unsigned RootFinder::MAXIT=50;

float
RootFinder::root(float x1, float x2)
{	// this code is almost exactly like rtsafe
	unsigned j;
	float df,dx,dxold,f,fh,fl;
	float temp,xh,xl,rts;

	funcd(x1,&fl,&df);
	funcd(x2,&fh,&df);
	if ((fl> 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
	   return FLT_MAX;
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if( fl < 0.0 ) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	
	}
	rts=(float)0.5*(x1+x2);
	dxold=(float)fabs(x2-x1);
	dx=dxold;
	funcd(rts,&f,&df);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=(float)0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		funcd(rts,&f,&df);
		if (f < 0.0) {
			xl=rts;
			fl=f;
		} else {
			xh=rts;
			fh=f;
		}
	}
	return FLT_MAX;

}


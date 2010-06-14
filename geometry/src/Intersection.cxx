// $Id$
//

#include "geometry/Intersection.h"
#include <cmath>
#include <cfloat>
#include <cassert>

static double step_parameter = 0.5;  // amount of curvature to step, in radians
static double tolerance = 2e-8;	    // how close is close enough (should be same as in Volume)

double
Intersection::distance(double maxStep, int sign)
{
#if 0 // play with this later
    double side = surf.how_near(ray.position()); // initial side
    assert( side*sign >= -1e-6);   // must be on correct side??
#endif

    // Step along the ray to find an intersection.
    // If sign + the leaving solutions are returned; if sign - entering sol. returned.

    double step = surf.distance(ray.position(), ray.direction(0), sign);

    // In case the ray is straight by-pass iterations and use exact solutions

    double curvature = ray.curvature();
    if(curvature == 0 ) return step ;

    // curved ray: make finite steps according to curvature until straight-line estimate is
    // within step
    double delta = step_parameter/fabs(curvature); // maximum step is 0.5 radian
    double start = 0;  // start of current step

    while( step > delta ) {
	if( (start += delta) > maxStep) {
	    // step is too far: if cut back, must still be on same side
	    double test = surf.how_near(ray.position(maxStep));
	    if( test*sign < 0 ) {
		// here if crossed boundary at maxStep, but distance missed
		start =	maxStep;
		step = surf.distance(ray.position(start), ray.direction(start), sign);
		break;
	    }
	    // conclude that surface is unreachable
	    return FLT_MAX;
	}
	step = surf.distance(ray.position(start), ray.direction(start), sign);
    }
    if(step==FLT_MAX) return FLT_MAX;

    // close: now iterate to precise solution
    while (fabs(step) > tolerance ) {
	start += step;
	double next_step = surf.distance(ray.position(start), ray.direction(start), sign);
	if( next_step == FLT_MAX) {
	    double test = surf.how_near(ray.position(maxStep));
	    assert(  test*sign > -tolerance);
	    return FLT_MAX; //give up
	}else if(next_step >0 && start> maxStep ) {
	    // prevent runaway
	    double test = surf.how_near(ray.position(maxStep));
	    assert(  test*sign > 0 );
	    return FLT_MAX;
	}

	if ( fabs(next_step) > fabs(step) )
	    next_step = step;
	step = next_step;
    }
    start += step;
#ifdef _DEBUG
    double test = surf.how_near(ray.position(start));
    assert( fabs(test) < 1e-6);
#endif
    return start;
}

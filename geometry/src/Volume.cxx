// $Id$
//
//

#include "geometry/Volume.h"
#include "geometry/Ray.h"
#include "geometry/Intersection.h"
#include <cassert>
# include <strstream>

static unsigned  s_count=0;
//THB 13-Aug-96 const double Surface_EPSILON= 2e-6; //(100*DBL_EPSILON);
const double Volume::Surface_EPSILON= 2e-5; //(100*DBL_EPSILON);

Volume::Volume(unsigned n)
: SurfaceList(n)
, max_dimension(FLT_MAX)
{
    s_count=0;
}

void Volume::addSurface(Surface* s)
{
    if( s_count<size() )
	operator[](s_count++)=s;
    else {
	push_back(s);
	s_count++;
    }
}

Volume::~Volume()
{  //  destructor: kill the surfaces
    deleteSurfaces();
}

void Volume::deleteSurfaces()
{  //  Removes all Surfaces from its list and deletes them
    for( iterator it=begin(); it !=end(); ++it) {
	delete *it; *it = 0;
    }
    s_count=0;
}

GeomObject&
Volume::transform(const CoordTransform& T)
{
   Shape::transform(T);
   for(unsigned i=0; i<surfaceCount(); i++)
	 surface(i).transform(T);
   return *this;
}

void Volume::printOn( std::ostream& os ) const
{
    Shape::printOn(os);
    os  << " maximum dimension   " << max_dimension   << "\n"
	<< "bounded by the surfaces:\n";
    for(unsigned i=0; i<surfaceCount(); i++)
	os << surface(i);
}

static unsigned isurf;   // local boundary index
int
Volume::getBoundaryIndex(const Point& p)const
{
    double closest = FLT_MAX;
    int found = -1;
    for( const_iterator it=begin(); it!=end(); ++it) {
	double d = fabs( (*it)->how_near(p) );
	if( d < closest ) {
	    found = it-begin();
	    closest = d;
	}
    }
    if( closest< Surface_EPSILON)
	return found;
    WARNING("Volume::associate: suface to tag not found");
    return -1;
}
int
Volume::lastBoundaryIndex()const
{
    // distanceToLeave sets index of left boundary: client can use

    return isurf;
}
// status of Ray::position() rel. to Surface
enum Surface_Status {Surface_INSIDE, Surface_AMBIGOUS, Surface_OUTSIDE} ;
/*
typedef int Surface_Status;
#define Surface_INSIDE      1
#define Surface_AMBIGOUS    2
#define Surface_OUTSIDE     3
*/
// type of solution of Intersection of Ray with Surface
enum Surface_Solution {Surface_SOLUTION, Surface_UNREACHABLE, Surface_NO_SOLUTION} ;
/*
typedef int Surface_Solution;
#define Surface_SOLUTION    1
#define Surface_UNREACHABLE 2
#define Surface_NO_SOLUTION 3
*/

static double Surface_distanceToLeave( const Surface &s, const Ray& r, double maxStep, Surface_Status& status, Surface_Solution& solution )
{
    double d = s.how_near( r.position() );

    if ( d > maxStep ) {
	// no solution due to max. step
		
	status   = Surface_INSIDE;
	solution = Surface_UNREACHABLE;
	return maxStep; //FLT_MAX/2;
    }


    if ( d < -Volume::Surface_EPSILON ) {
	// the outside case for s
		std::ostrstream error_text;
		error_text << "Volume--bad geometry for surface \n\t" << s
		    << "\n\t" << r.position() << " is " << d << "cm outside\n";
		error_text << '\0';
		WARNING(error_text.str());
		status   = Surface_OUTSIDE;
		solution = Surface_NO_SOLUTION;
		return FLT_MAX;
    }
    else if (d < Volume::Surface_EPSILON ) {
		// the ambigous case for s --> take first leaving intersection
		
	status   = Surface_AMBIGOUS;
	solution = Surface_SOLUTION;

#if 0 // THB: this screws me up
	if (s.normal( r.position() ) * r.direction() > 0)
	    return 0;
#endif

	Intersection inter(r,s);
	d = inter.distance(maxStep, 1);

	solution = (d == FLT_MAX) ? Surface_NO_SOLUTION : Surface_SOLUTION;
	return d;
    }
    else {
	// the inside case for s --> take first leaving intersection

	status   = Surface_INSIDE;

	Intersection inter(r,s);
	d = inter.distance(maxStep, 1);

	solution = (d == FLT_MAX) ? Surface_NO_SOLUTION : Surface_SOLUTION;
	return d;
    }
}


double
Volume::distanceToLeave( const Ray& r,  double maxStep ) const
{
	//  Distance along the Ray r or any of its derived classes
	//  to the nearest Surface to leave the Volume, constrained to be
	//  less than maxStep.
	//  If no intersection point is found, returns a large number (FLT_MAX).
	//  Loop over each Surface in the Volume's list.
	//  Use the how_near function to find the closest distance between the
	//  origin of the Ray (or its derived class) and the Surface.
	//  Only if that distance is non-negative, less than maxStep,
	//  and less than the distance along the Ray found for the previous Surface,
	//  find the distance along the Ray (or its derived class) to that Surface.
	//  This will get the shortest distance along the Ray to leave the Volume.
	//  The distanceToLeaveSurface member function of Ray will make sure
	//  that the Ray or its derived class is leaving the Surface (defined by
	//  the direction of the normal to the Surface) at the point of intersection.

    double             t = FLT_MAX;
    Surface_Status     statusVolume   = Surface_INSIDE;
    Surface_Solution   solutionVolume = Surface_NO_SOLUTION;

    isurf = 99999;
    for (unsigned i=0; i< surfaceCount(); i++) 	{
	const Surface& s = surface(i);
	Surface_Status   status;
	Surface_Solution solution;
	double d = Surface_distanceToLeave( s, r, maxStep, status, solution );

	if ( status == Surface_OUTSIDE )
	{
	    statusVolume   = Surface_OUTSIDE;
	    solutionVolume = Surface_NO_SOLUTION;
	    return FLT_MAX;
	}

	// special speedup: when one surface finds the current point to be on-surface and the ray
	// leaving that surface at the current point, then there is no need to test all surfaces
	// again for the inside condition, because the current point has been tested against all
	// other surfaces already and is really inside or on-surface
	if (t == 0)
	{
	    statusVolume   = Surface_AMBIGOUS;
	    solutionVolume = Surface_SOLUTION;
	    return 0;
	}

	if ( status > statusVolume )
		statusVolume = status;

	if ( solution < solutionVolume )
		solutionVolume = solution;

	if (d < t)
	{
	    isurf = i;
	    t = d;
	}
    }


    // at this point, i know that i'm INSIDE or ON-SURFACE relativ to ALL surfaces !!!
    // It might be, that there was no intersection found, though. there may be 3 reasons:
    //   (a) the surface which i would have hit was unreachable
    //   (b) the ray is circling
    //   (c) Volume is has no closed surface (e.g. a box with one side missing)
    //       this is clearly an error
    if (solutionVolume == Surface_UNREACHABLE)
	// (a) there still might be an error in the geometry of Volume, but
	//     there is no easy way to find out !!
	return t; //FLT_MAX/2;

   if (solutionVolume == Surface_NO_SOLUTION) {
	// cases (b) and (c), but how am i supposed to find out, which one ?
#if 0 // do not need anymore? def _DEBUG
       std::cerr << "Volume::distanceToLeave, can't leave!\n";
       std::cerr << "\tVolume=" << *this << '\n';
       std::cerr << "\tRay= " << r << '\n';
#endif
	return FLT_MAX;
   }

    // maybe another speedup (?):
    // we already know we're inside. so, the closest point p where the ray leaves MUST be
    // inside (or rather on-surface), because if it would be outside the ray must have
    // left the volume earlier (if we have a contigous ray that is), so we should have
    // found a smaller solution in contradiction to the assumption of p being the closest
    // solution.
    //
    // PROBLEM:  what about the right surface, which would give the real, closest
    // point being unreachable?
    //	there are examples, where the surface
    //
    // so maybe it's safe to return here. NO !
    // return t;


    // check for r.position(t) to be inside or on-surface
    Point p = r.position(t);
    int bad=0;
    for (unsigned j=0; j<surfaceCount(); j++)
    {	
	if (j==isurf)
		continue;

	bad = (surface(j).how_near(p) < -Surface_EPSILON);

	if ( bad )
	{
	    std::ostrstream cerr;
	    cerr << "Volume::distanceToLeave, position bad by "
		<< surface(j).how_near(p) << '\n';
	    cerr << "\tVolume=" << *this << '\n';
	    cerr << "\tRay= " << r << '\n';
	    cerr << "\tPoint" << p << '\n';
	    cerr << '\0'; WARNING(cerr.str());
	    solutionVolume = Surface_NO_SOLUTION;
	    return FLT_MAX;
	}
    }

    return t;
}


double
Volume::distanceToEnter( const Ray& r, double maxStep ) const
{ //  Distance along the Ray r or any of its derived classes
  //  to the nearest Surface of the Volume, constrained to be < maxStep.
  //  If the Ray's origin is inside and no entering solution found: return -1.
  //  If no intersection point is found: return a large number (FLT_MAX).
  //  Quick test to see if the Ray can possibly enter this Volume:
  //  if the distance from the origin of the Ray to the center of the Volume
  //  minus the maximum dimension (which is the radius of a sphere which
  //  circumscribes the Volume) is greater than maxStep, then the Ray cannot
  //  enter the Volume.

    Point startPos = r.position();

    if( maxStep < FLT_MAX ) {
	Vector oc = center() - startPos;
	if ( ( oc.mag() - max_dimension ) > maxStep ) return FLT_MAX;
    }

    double t = FLT_MAX;             // value to be returned

    for(unsigned i=0; i<surfaceCount(); i++) {

        //  Check that we can reach surface in this step: examine signed
	//  closeness of the surface
	double d = surface(i).how_near(startPos);
	if( d> Surface_EPSILON	    // already inside: can't enter thru this surface
	    ||  fabs(d) > t	    // alredy have a closer surface
	    ||  fabs(d) > maxStep   // no point in trying
	 )	continue;

    // Now check distance to the surface: create an Intersection object
	Intersection inter(r,surface(i));
	d = inter.distance(maxStep, -1 );
	if( d< 0 && d > -Surface_EPSILON)
		d=0;  // right on this surface!!
		
	if( d < t && d >= 0 ){

    // It's the closest, so far, but it must be inside the others
	    Point p = r.position(d);
	    int bad=0;
	    for(unsigned j=0; j<surfaceCount(); j++) {
		if(j==i) continue;
		if( (bad = !surface(j).inside(p))!=0)break;
		    }
		if( !bad ) {
		    t = d; // must be the one.
			break;
	    }
        }
    }
    return t;
}

int
Volume::inside( const Point& x ) const
{
    //  Return 0 if a point x is outside the Volume, 1 if inside.
    //  The point is inside the Volume if it is inside all the Surfaces
    //  which define the Volume.
    //  Surfaces provide an outward-pointing unit normal for this purpose.
    for(unsigned i=0; i<surfaceCount(); i++)
    {
	if( !surface(i).inside( x ))
	   return 0;
    }
    return 1;
}

double Volume::howNear( const Point& x) const
{
   double r=FLT_MAX;
    for(unsigned i=0; i<surfaceCount(); i++) {
	double d =surface(i).how_near( x );
        if (d < r) r=d;
    }
    return r;
}
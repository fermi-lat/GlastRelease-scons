// $Id$
//
//

#ifndef GEOMETRY_SURFACE_H
#define GEOMETRY_SURFACE_H

#include "geometry/Point.h"
#include <float.h> // for FLT_MAX on NT
/**
 Virtual class that defines properties of an abstract surface
 Contains two data members:
<UL>
<LI>a reference to a Point corresponding to an origin that it does not own,
    and may be sharing with other Surfaces
<LI> a Vector reference direction, which may be interpreted by subclasses as desired
</UL>
*/
class Surface : public GeomObject
{

 public:
 
    /// return true if point is "inside" the directed surface default is just the inline above
    virtual int inside(const Point& p)const{ return ( how_near( p ) >= 0 ); }


    /** signed function giving minimum distance to the surface: surface
     is actually defined by zero.
     The distance will be positive if the point is inside the Plane,
     negative if the point is outside.
    */
    virtual double how_near( const Point& x ) const=0;

 
    /// distance along direction v from point x to enter or leave the surface.
    virtual double distance( const Point& x, const Vector& v, int inout)const=0;

    /// returns normal-vector at Point x, assuming x to be on the surface
    virtual Vector normal( const Point& x )const=0;

    /// apply transformation to referenceDirection. 
    /**Subclass must apply to other directions
     note that a Surface will not respond to a translation transformation.
    */
    virtual GeomObject& transform(const CoordTransform& T );

    /// reverse sense of surface
    virtual void reverse();

    /// access to origin
    const Point& origin()const{return referencePoint;}

    /// access to direction
    const Vector& direction() const { return referenceDirection; }

  
    /// virtual copy, must be implemented by subclass to allow copy
    virtual Surface& copy()const;

    /** must be supplied by subclass if appropriate. Allows the surface to be "moved" by
     changing the parameter describing the offset of the surface from the origin. 
    */
    virtual void setOffset(double){};
    virtual double offset()const{return FLT_MAX;};


 protected:
    Surface( const Point& o, const Vector& n )
	: referencePoint( o ), referenceDirection(n) 
    {}
    // only subclasses can implement

    Surface(const Surface&);
    // copy constructor not accessible

 
 private:

    const Point& referencePoint;  
    // reference point to use: defined elsewhere
  
    Vector referenceDirection;    
    // direction: normal for a plane surface, axis  for cylinder


};
//---------------------------------------------------------------------
//                 inlines
//---------------------------------------------------------------------

inline std::ostream& operator<<( std::ostream& os, const Surface& s )
  { s.printOn(os); return os;}

inline std::ostream& operator<<( std::ostream& os, const Surface* s )
  { s->printOn(os); return os;}

#endif


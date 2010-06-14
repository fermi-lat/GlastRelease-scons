// $Heading: Track.h $
//   Author: T. Burnett
//

#ifndef GEOMETRY_TRACK_H
#define GEOMETRY_TRACK_H

#include "geometry/GeomObject.h"
#include "geometry/Point.h"
class Ray;
#include <string> // gcc bug
#include <vector>

/// A Track is a list of Ray objects: it will document the complete trajectory of a particle.
class Track  : public GeomObject
{
public:
    /// constructor: must specify first ray and if not charged
    Track( Ray* first, bool charged=true);
    /// destructor: deletes all rays in the list
    ~Track();

    GeomObject& transform(const CoordTransform&);

    /// add a new ray: assumes ownership
    void addSegment( Ray* next );

    /// position and direction
    Point position( double s ) const;
    Vector direction( double s ) const;

    bool charged()const{return m_charged;}

    const char *nameOf() const { return "Track"; }
    void printOn( std::ostream& os = std::cout ) const;


  typedef std::vector<Ray*> Raylist;

    typedef Raylist::const_iterator const_iterator;
    typedef Raylist::iterator iterator;
    const_iterator begin()const{return rayList.begin();}
    const_iterator end()const{return rayList.end();}
    iterator begin() {return rayList.begin();}
    iterator end() {return rayList.end();}

private:
   void findSegment( double s );
   // used to find the right segment

   static Ray* currentRay;
   static double currentArcLength;
   // set by findSegment


   static float maxKink;
   // maximum angle tolerated in line segment approximation to
   // curved tracks.

    bool m_charged;	// true if charged
   float arclength;	 // length along the Track


   Raylist rayList; // the list of Rays
};

#endif


// $Heading: Track.h $
//   Author: T. Burnett
//
// A Track is a list of Ray objects: it will document the complete trajectory of
// a particle.

#ifndef __TRACK_H
#define __TRACK_H

#include "geometry/GeomObject.h"
#include "geometry/Point.h"
class Ray;
#include <string> // gcc bug
#include <vector>


class Track  : public GeomObject
{
public:
    Track( Ray* first, bool charged=true);
    // constructor: must specify first ray and if not charged
    ~Track();
    // destructor: deletes all rays in the list

    GeomObject& transform(const CoordTransform&);

    void addSegment( Ray* next );
    // add a new ray: assumes ownership

    Point position( double s ) const;
    Vector direction( double s ) const;
    // position and direction

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


// $Heading: Track.cxx $
//   Author: T. Burnett
//

#include "geometry/Track.h"
#include "geometry/Ray.h"


Ray* Track::currentRay;
double Track::currentArcLength;


float Track::maxKink = 0.2f;//0.05;

Track::Track( Ray* first, bool charged)
: m_charged(charged)
{
   rayList.push_back( first);
   arclength = first->getArcLength();
   currentRay = first;
}

Track::~Track()
{
    Raylist::iterator rayptr = rayList.begin();
    while( rayptr < rayList.end() )
	delete *rayptr++;
}

GeomObject&
Track::transform(const CoordTransform& T)
{
  for(unsigned i=0; i< rayList.size(); i++)
    rayList[i]->transform(T);
  return *this;
}
void Track::findSegment( double s )
{  unsigned i=0;
   currentRay = rayList[0];
   currentArcLength = 0.;
   double ds = currentRay->getArcLength();
   while ( i<rayList.size() && s > currentArcLength+ds) {
	currentArcLength += ds;
	currentRay =rayList[++i];
	ds = currentRay->getArcLength();
   }
}

void Track::addSegment( Ray* next )
{
  float ds = next->getArcLength();
  arclength += ds;
/*
//  This section of code fails when multiple scattering is present
   if( currentRay->curvature()==0
	    && next->curvature()==0 )
   {
      // this is a continuation of a straight ray. Eventually want
      // to generalize it to helices
       currentRay->addArcLength( ds );
       delete next;
   }
   else
   {   append( (GmsListLink*)next );
       currentRay = next;
   }
*/
//   Till I(WBA) have time to fix this
   rayList.push_back( next );
   currentRay = next;

}

void Track::printOn(std::ostream& os) const
{
    os << "Track:  total arclength = "<< arclength << '\n';
    for( unsigned i=0; i<rayList.size(); i++ )
	rayList[i]->printOn(os);
}

Point
Track::position( double s ) const
{
   Track* self = (Track*)this; // have to cast to non-const
   self->findSegment(s);
   return currentRay->position(s-currentArcLength);
}

Vector Track::direction( double s ) const
{
   Track* self = const_cast<Track*>(this); // have to cast to non-const
   self->findSegment(s);
   return currentRay->direction(s-currentArcLength);
}




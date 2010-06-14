//  $Heading: PointRep.h $
//   Author: T. Burnett

#ifndef POINTREP_H
#define POINTREP_H

#include "graphics/GraphicsRep.h"

class Hep3Vector;

class PointRep : public GraphicsRep 
// class that represents a single point in ArveGraphics
{

public:
    PointRep(const Hep3Vector& t);
    // constructor: save reference to given point

    virtual void update();
    // create a representation of the point

private:
    const Hep3Vector & m_point;
};
// Inlines
inline PointRep::PointRep(const Hep3Vector& t):m_point(t){}

#endif // POINTREP_H



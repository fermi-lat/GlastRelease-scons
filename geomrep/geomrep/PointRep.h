//  $Heading: PointRep.h $
//   Author: T. Burnett

#ifndef POINTREP_H
#define POINTREP_H

#include "gui/DisplayRep.h"

namespace CLHEP {class Hep3Vector;}

class PointRep : public gui::DisplayRep 
// class that represents a single point in ArveGraphics
{

public:
    PointRep(const CLHEP::Hep3Vector& t);
    // constructor: save reference to given point

    virtual void update();
    // create a representation of the point

private:
    const CLHEP::Hep3Vector & m_point;
};
// Inlines
inline PointRep::PointRep(const CLHEP::Hep3Vector& t):m_point(t){}

#endif // POINTREP_H



//     $Id$
//  Author: T. Burnett
// Project: Arve graphics 

#ifndef ARCREP_H
#define ARCREP_H
#include "graphics/ArvePlottableRep.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"

class ArcRep : public ArvePlottableRep {
public:
    ArcRep( const Point& center, double radius, const Vector& xhat,
	const Vector& yhat, double phi0);
    // constructor: specify parameters

    void update();
    
private:
    Point m_center;
    double m_radius;
    Vector m_xhat, m_yhat;
    double m_phi0;
};

#endif //ARCREP_H

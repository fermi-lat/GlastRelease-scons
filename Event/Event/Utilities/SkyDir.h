// $Header$
#ifndef EVENT_SKYDIR_H
#define EVENT_SKYDIR_H 1


// Include files
#include <iostream>
#include "geometry/CoordTransform.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

/*! \class SkyDir
\brief encapsulate a direction.
This holds (and returns) an absolute direction, in terms of inertial coordinate systems
*/
enum coordSystem { 
        GALACTIC,  //!  fixed direction with respect to the galactic coordinate system (l,b)
        CELESTIAL //! fixed direction with respect to the celestial coordinate system (ra,dec) in the J2000 epoch.
};

class SkyDir
{
public:
    ///Constructors
    //(l,b) or (Ra, Dec) instantiation
    SkyDir(double param1, double param2, coordSystem inputType = GALACTIC);
    SkyDir(Hep3Vector);
    
    ///return methods
    Hep3Vector& operator () () {return m_dir;}
    double l () const;
    double b () const;
    double ra () const;
    double dec () const;
    Hep3Vector r () const;
    Hep3Vector m_dir;

private:
    Rotation celToGal() const;
    std::pair<double,double> setCelCoordsFromDir() const;
    std::pair<double,double> setGalCoordsFromDir() const;

};

#endif    // LHCBEVENT_SKYDIR_H

// $Header$
#ifndef EVENT_SKYDIR_H
#define EVENT_SKYDIR_H 1


// Include files
#include <iostream>
//#include "GaudiKernel/StreamBuffer.h"
//#include "Event/TopLevel/Definitions.h"
#include "geometry/Vector.h"
#include "geometry/CoordTransform.h"
#include "CLHEP/Vector/Rotation.h"

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
    SkyDir(Vector);
    
    ///return methods
    double l();
    double b();
    double ra();
    double dec();
    Vector r();
    

private:
    Rotation celToGal();
    void setCelCoordsFromDir();
    void setGalCoordsFromDir();
    double m_l,m_b,m_ra,m_dec;
    Vector m_dir;

};


  /*  {

public:

  /// Constructors
  TimeStamp()
    : m_time(0)                                                              { }
  TimeStamp( double t )
    : m_time(t)                                                              { }
  /// Destructor
  ~TimeStamp()                                                               { }

  /// Retrieve time
  double time() const                                                            {
    return m_time;
  }
  /// Update time 
  void setTime( double value )                                                   {
    m_time = value;
  }

  operator double()const { return time(); }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s, const TimeStamp& obj )    {
    return s << obj.m_time;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s, TimeStamp& obj )          {
    return s >> obj.m_time;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const TimeStamp& obj )    {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
  std::ostream& fillStream( std::ostream& s ) const                            {
    return s << "class TimeStamp : "
      << EventField( EventFormat::field12 )
      << m_time;
  }

private:

  /// Time
  double m_time;

};
*/

#endif    // LHCBEVENT_SKYDIR_H

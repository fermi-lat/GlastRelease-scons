// $Id$
/** 
* \class StaticOrbit
*
* \brief Child-class of the Orbit class which defines a static - specifiable
* orbit designed for situations in which the orbit is specified by some
* external data source.
*
* \author Sawyer Gillespie - hgillesp@u.washington.edu
* 
* $Header $
*/


#include "flux/Orbit.h"


class StaticOrbit : public Orbit {
public:
    StaticOrbit (double lt = 0., double ln = 0.) : m_lat(lt), m_lon(ln) {}
    
    // overloads
    virtual double  latitude ( double time ) { return m_lat; }
    virtual double  longitude ( double time ) { return m_lon; }
    virtual double  pitch ( double time ) { return m_pitch; }
    virtual double  yaw ( double time ) { return m_yaw; }
    virtual double  roll ( double time ) { return m_roll; }
    
protected:
    // set methods
    virtual void  setLatitude ( double l ) { m_lat = l; }
    virtual void  setLongitude ( double l ) { m_lon = l; }
    virtual void  setPitch ( double p ) { m_pitch = p; }
    virtual void  setYaw ( double y ) { m_yaw = y; }
    virtual void  setRoll ( double r ) { m_roll = r; }
    
private:
    double m_lat;   // latitude of the instrument
    double m_lon;   // longitude of the instrument
    double m_pitch; // pitch of the instrument
    double m_yaw;   // yaw of the instrument
    double m_roll;  // roll of the instrument
};
// $Id$


#ifndef ORBIT_H
#define ORBIT_H
/** 
* \class Orbit
*
* \brief Calculates the position of a satellite in a low, circular orbit.
*
* The path is assumed to go to the east, with the rotation of the earth,
* as is the case for most satellite orbits.
* All angles are expressed in degrees.  Longitude increases to the east.
* <hr>
* Constructor:
* The three arguments are
* 1.  Longitude of the ascending node (where the orbit crosses
* the equator, moving northward) in degrees
* 2.  Altitude in kilometers.  Defaults to 600.
* 3.  Orbital inclination in degrees.  Defaults to 28.5.
* <hr>
* Member functions:
* latitude(time) and longitude(time) give the position of the
* satellite (in degrees) when "time" minutes have elapsed since 
* the ascending node passage.  longitude is between 0 and 360.
* coords(time) returns the latitude and longitude as a pair<double,double>.
* period() returns the period of the orbit (in minutes).
* inclination(), altitude(), and ascendingLon() return the parameters
* used in the constructor, in case you forget.
* There is no need for an explicit destructor, copy constructor,
* or assignment operator.
* 
* $Header $
*/

#include <utility>        // for STL pairs
#include "geometry/CoordTransform.h"

class Orbit{
    
public:
    /// default constructor
    Orbit(double asclon = 0., double alt = 600., double inc = 28.5, double phs = 0.,
        double apitch = 0., double ayaw = 0., double aroll = 0.);
    
    /// phase of the orbit as a proportion of a single period
    double get_phase (double time) const;
    
    ///returns position as a pair using  specified time
    std::pair<double,double> coords(double time) const;

    /*! 
    \param latitude Latitude is in the range (-180,180)
    */
    /// latitude as a function of time (in minutes)    
    virtual double latitude(double time) const;
    
    /*! 
    \param longitude Longitude is in the range (0,360)
    */
    /// longitude as a function of time, taking into account the starting
    /// longitude and the eastward rotation of the earth
    virtual double longitude(double time) const;
    
    /// pitch of the spacecraft (rotation around E-W axis)
    /// zenith pointing is pitch == 0
    virtual double pitch (double time) const;
    
    /// yaw of the spacecraft (rotation around zenith)
    /// square to E-W and N-S axis represents yaw == 0
    virtual double yaw (double time) const;
    
    /// roll of the spacecraft (rotation around N-S axis)
    /// zenith pointing is roll == 0
    virtual double roll (double time) const;
    
    /// phase of the orbit at a given time
    virtual double phase (double time) const;
    
    /// set the inclination of the orbit
    void    inclination ( double inc );
    
    /// set the longitude of the ascending node.
    void    ascendingLon ( double asc );
    
    /// set the phase of the orbit (in radians)
    void    startphase ( double phi );
    
    /// explicitly output the rotation matrix to the screen
    void displayRotation(Rotation rot);
    
    /// inlined access methods
    inline double startphase() const;
    inline double period() const;
    inline double inclination() const;
    inline double altitude() const;
    inline double ascendingLon() const;
    
    //transformations for the current position (to the zenith-pointing coordinate system)
    Rotation CELTransform(double time);
    
    
    Rotation Orbit::latLonTransform(double time) const;
    
    double Orbit::testLongitude(double time) const;
    double Orbit::testLatitude(double time) const; 
    
protected:
    // manipulation of the orbital parameters - defaults do nothing, however
    // subclasses may overload depending upon which orbital parameters they 
    // specify
    
    virtual void setLatitude ( double ) {}
    virtual void setLongitude ( double ) {}
    virtual void setPitch ( double ) {}
    virtual void setYaw ( double ) {}
    virtual void setRoll ( double ) {}
    
    
    // computation of pointing characteristics of GLAST - assumed, here, to be zenith-pointing
    void computeAttitudes(double time);
    
    // friends - so that only the boss can manipulate this class
    friend class GPS;
    
private:
    double m_ascendingLon;   // longitude of ascending node (degrees)
    double m_inclination;    // orbit inclination (degrees)
    double m_altitude;       // altitude (km)
    double m_period;         // orbital period (minutes)
    double m_sini;           // sine of inclination angle
    double m_cosi;           // cosine of inclination angle
    double m_startphase;     // starting phase offset of the orbit
    double m_precessPeriod;   // Period of orbital precession
    
    double m_degsPerRad;
    double m_secsperday;
    
    
    //RA, DEC of GLAST X and Y axes:
    double m_rax;
    double m_raz;
    double m_decx;
    double m_decz;
    
};

inline double Orbit::period() const {return m_period;}

inline double Orbit::inclination() const {return m_inclination;}

inline double Orbit::altitude() const {return m_altitude;}

inline double Orbit::ascendingLon() const {return m_ascendingLon;}

inline double Orbit::startphase () const {return m_startphase;}

#endif ORBIT_H

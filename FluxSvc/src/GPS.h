// $Header$

#if !defined(_H_GPS_CLASS)
#define _H_GPS_CLASS

/** 
* \class GPS
*
* \brief Models the Global Positoning System for a spacecraft. Handles time, position, and orientation for the instrument as a whole.
* 
* $Header $
*/


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "facilities/Scheduler.h"
#include "facilities/Observer.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "geometry/CoordTransform.h"

#include <iostream>

class Orbit;
// 
//  class GPS
//
//!  Represents the Global Positioning System on-board the spacecraft. An Orbit
/*!  object is used to compute the spacecraft's position and pointing characteristics.
Time is tracked through this object, and synchronized with the Scheduler for 
discrete event simulation. An expansion factor is provided to allow for acceleration
or randomization of the orbit. If the expansion factor is negative, then the position
of the spacecraft is chosen as a random distribution over an orbit. Otherwise, the expansion
factor represents an acceleration of the spacecraft's orbit. Ie. an expansion factor
of 2 would reduce the orbit period of the spacecraft by 1/2.
*/
class GPS  
{
public:
    class Coords {
    public:
        Coords( double alat, double alon, double apitch
            , double ayaw, double aroll, GPStime atime, double phase ); 
        Coords();
        
        GPStime time () const { return m_time; }
        double  lat () const { return m_lat; }
        double  lon () const { return m_lon; }
        double  pitch () const { return m_pitch; }
        double  yaw () const { return m_yaw; }
        double  roll () const { return m_roll; }
        double  phase () const { return m_phase; }
        
    private:
        GPStime m_time;
        double m_lat, m_lon, m_pitch, m_yaw, m_roll, m_phase;
    };
    
    // const access
    
    /// GPS synchronized time for the satellite
    GPStime	time () const; 
    /// present latitude
    double	lat () const; 
    /// present longitude
    double	lon () const; 
    /// present 'pitch' angle in N-S direction (roll around x-axis)
    double	pitch () const; 
    /// present 'yaw' angle - (rotate around z-axis)
    double	yaw () const; 
    /// present 'roll' andle in E-W direction (roll around y-axis)
    double      roll () const; 
    /// get the current phase of the orbit
    double      phase () const;
    /// access to the orbit (for computation only)
    const Orbit* orbit () const;
    /// current state of the GPS system
    Coords      state () const;  
    /// expansion of the current orbit
    double      expansion () const; 
    /// sample interval for random orbit distribution
    GPStime     sampleintvl () const; 
    /// return the orbit's ascending longitude 
    double     ascendingLon()const; 
    /// access m_rotangles
    std::pair<double,double> rotateAngles(); 
    // set data
    
    /// set a specific Orbit object for lat/lon calculations
    void    orbit ( Orbit* );
    /// pass a specific amount of time
    void    pass ( double );
    /// set the current phase of the orbit
    void    phase ( double );
    /// set the expansion factor for the orbit (-1) = random
    void    expansion ( double );
    /// synchronize w. scheduler
    void    synch ();    
    /// set the sample interval
    void    sampleintvl ( GPStime );
    /// special to set the ascending longitude	
    void    ascendingLon(double);   
    /// set m_rotangles
    void    rotateAngles(std::pair<double,double> coords); 
    
    /// print time & position
    void    printOn(std::ostream& out) const; 
    
    // notification
    void notifyObservers() { m_notification.notify();}
    Subject&    notification() { return m_notification; }
    
    // static access/destruction
    static GPS*	instance();
    static void     kill ();
    
    ///return galactic-pointing coordinates of the GLAST satellite
    std::pair<double,double> GPS::galPositionOfGlast();
    
    /// transform functions using galPositionOfGlast.  This infrastructure is no longer used.
    Vector earthToGlast(Vector launchDir);
    Vector galaxyToGlast(Vector launchDir);
    
    /// return the rotation for compensation for the rocking angles.
    Rotation rockingAngleTransform(double time);
    
    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    Rotation GPS::transformGlastToGalactic(double time);
    
    
    
    Orbit*  orbit ();               // access the orbit (for manipulation)
    protected:
        // singleton - protect ctor/dtor
        GPS();
        virtual ~GPS();
        
        // only class/subclasses have access
        void    pitch ( double );       // set pitch
        void    yaw ( double );         // set yaw
        void    roll ( double );        // set roll
        void    lat ( double );         // set latitude
        void    lon ( double );         // set longitude
        void    time ( GPStime );       // set time
        GPStime orbittime () const;     // access to orbit time
        void    orbittime ( GPStime );  // set orbit time
        void    setState ( const GPS::Coords& ); // set the orbital parameters
        std::pair<double,double> m_rotangles;  //angles for coordinate rotation (rocking angle)
        
        
        
        
        // friends
        friend class FluxGenerator;
        
    private:
        static GPS* s_instance;
        
        Orbit*  m_orbit;        // orbital position object 
        double  m_expansion;    // orbit expansion factor
        GPStime m_time;	    // global time
        GPStime m_orbittime;    // time in orbit (depends on expansion)
        double  m_sampleintvl;  // interval to sample for each pt. in the orbit - to normalize spectra
        
        // notification
        Subject    m_notification; 
        
};

inline std::istream&    operator>>(std::istream& i, GPS::Coords& c) {
    double  p, y, r, lt, ln, tm, ph;
    i >> lt; i >> ln; i >> p; i >> y; i >> r; i >> tm; i >> ph;
    c = GPS::Coords( lt, ln, p, y, r, tm, ph );
    return i;
}

inline std::ostream&    operator<<(std::ostream& o, const GPS::Coords& c) {
    o << ' ' << c.lat() << ' ' << c.lon() << ' ' 
        << c.pitch() << ' ' << c.yaw() << ' ' 
        << c.roll() << ' ' << c.time() <<' ' << c.phase();
    return o;
}
#endif // !defined(AFX_GPS_H__F9844433_4E64_11D2_B4DD_00A0C9960210__INCLUDED_)

// GPS.cxx: implementation of the GPS class.
// $Id$
//////////////////////////////////////////////////////////////////////

#include "GPS.h"

#include "Orbit.h"
#include "CLHEP/Random/RandFlat.h"

#include <iomanip>


// static declaration

GPS*	GPS::s_instance = 0;

GPS::GPS() 
: m_orbit(new Orbit),
m_rockDegrees(25.),
m_rockType(NONE),
m_earthOrbit(new astro::EarthOrbit),
m_expansion(1.),    // default expansion:regular orbit for now
m_time(0.), 
m_orbittime(m_time),
m_sampleintvl(0.001),
m_rotangles(std::make_pair<double,double>(0.,0.))
{}

GPS::Coords::Coords( double alat, double alon, double apitch
                    , double ayaw, double aroll, GPStime atime, double aphase) 
                    : m_time(atime), m_lat(alat), m_lon(alon), 
                    m_pitch(apitch), m_yaw(ayaw), m_roll(aroll), m_phase(aphase)
{}

GPS::Coords::Coords(){}


GPS::~GPS ()
{ delete m_orbit; }


double	GPS::pitch () const
{ 
    return orbit()->pitch( orbittime() ); 
}

double	GPS::yaw () const
{ 
    return orbit()->yaw( orbittime() );   
}

double   GPS::roll () const
{ 
    return orbit()->roll( orbittime() ); 
}

const Orbit*   GPS::orbit () const 
{
    return m_orbit;
}

Orbit*   GPS::orbit () 
{
    return m_orbit;
}

void GPS::synch ()
{
    static bool first=true;
    bool changed=  false;
    
    if (Scheduler::instance()->running()) {
        time( Scheduler::instance()->elapsed_time() );
        changed = true; // maybe have threshold before nofitying?
    }
    if (expansion() < 0.) {
        static GPStime  last_time = time();
        if ((time() - last_time) > m_sampleintvl) {
            orbittime( RandFlat::shoot(orbit()->period()) );
            orbit()->ascendingLon(RandFlat::shoot(360.));
            last_time = time();
            changed = true;
            
        }
    }
    else if( expansion()>0 ) {
        orbittime( time()*expansion() );
        changed=true; 
    }
    
    // notify observers if changed (or first time thru)
    if( changed || first) notifyObservers();
    first=false;
    
}

double   GPS::lat () const
{
    return orbit()->latitude( orbittime() ); 
}

double	GPS::lon () const
{ 
    return orbit()->longitude( orbittime() ); 
}

GPStime	GPS::time ()  const
{ 
    return m_time;
}

double   GPS::phase () const
{
    // compute orbital phase based upon current time
    return orbit()->phase( orbittime() );
}

double   GPS::expansion () const
{
    return m_expansion;
}

void	GPS::orbit ( Orbit* o )
{
    if (o == 0) return;
    delete m_orbit;
    m_orbit = o;
}


void GPS::pass ( double t )
{ 
    if (!Scheduler::instance()->running())	{
        time(time() + t);
    }   // cannot pass when scheduler is in control!
}

void GPS::phase ( double p )
{
    // jump to a particular point in the orbit
    double  pd = orbit()->period();
    if (pd != 0.) orbit()->startphase(p - 2.*M_PI*time()/pd);
}

void GPS::expansion ( double e )
{
    m_expansion = e; 
}

void GPS::pitch ( double p )
{
    m_orbit->setPitch(p);
}

void GPS::yaw ( double y )
{
    m_orbit->setYaw(y);
}

void GPS::roll ( double r )
{
    m_orbit->setRoll(r);
}

void GPS::lat ( double l )
{
    m_orbit->setLatitude(l);
}

void GPS::lon ( double l )
{
    m_orbit->setLongitude(l);
}

void GPS::time ( GPStime t )
{
    m_time = t;
}

GPS*	GPS::instance() 
{ return (s_instance != 0) ? s_instance : (s_instance = new GPS()); }

void GPS::kill()
{
    delete s_instance;
    s_instance = 0;
}

GPS::Coords GPS::state () const
{
    return GPS::Coords( lat(),lon(),pitch(),yaw(),roll(),time(), orbit()->phase(time()) );
}

void GPS::orbittime ( GPStime t ) 
{
    m_orbittime = t;
}

GPStime GPS::orbittime () const
{
    return m_orbittime;
}

void    GPS::sampleintvl ( GPStime s )
{
    m_sampleintvl = s;
}

GPStime  GPS::sampleintvl () const
{
    return m_sampleintvl;
}

void    GPS::setState ( const GPS::Coords& c )
{
    m_orbit->setLatitude(c.lat());
    m_orbit->setLongitude(c.lon());
    m_orbit->setPitch(c.pitch());
    m_orbit->setYaw(c.yaw());
    m_orbit->setRoll(c.roll());
    if (c.time() > 0.) time(c.time());  // if the time is relevant, use that time
}
void GPS::ascendingLon(double lon)
{
    m_orbit->ascendingLon(lon);
}
double GPS::ascendingLon()const
{
    return m_orbit->ascendingLon();
}
void    GPS::printOn(std::ostream& out) const
{
    out << "GPS: time=" << time() 
        << ", lat, lon=" 
        << std::setprecision(4) << lat() << ", " << lon() 
        << std::endl;
}


//access m_rotangles
std::pair<double,double> GPS::rotateAngles(){
    return m_rotangles;
    
}

//set m_rotangles
void GPS::rotateAngles(std::pair<double,double> coords){
    m_rotangles=coords;
}


Rotation GPS::rockingAngleTransform(double seconds){
    //Purpose:  return the rotation to correct for satellite rocking.
    //Input:  Current time
    //Output:  3x3 rocking-angle transformation matrix.
    using namespace astro;
    
    double time = m_earthOrbit->dateFromSeconds(seconds);
    
    double inclination = m_earthOrbit->inclination();
    double orbitPhase = m_earthOrbit->phase(time);
    m_position = m_earthOrbit->position(time);
    
    SkyDir dirZ(m_position.unit());
    SkyDir dirX(dirZ.ra()-90., 0.0);
    
    //rotate the x direction so that the x direction points along the orbital direction.
    dirX().rotate(dirZ.dir() , inclination*cos(orbitPhase));
    
    //now set the zenith direction to do the rocking properly.
    m_RAZenith = dirZ.ra();
    m_DECZenith = dirZ.dec();

    double rockNorth = m_rockDegrees*M_PI/180;
    
    //here's where we decide how much to rock about the x axis.  this depends on the 
    //rocking mode.
    if(m_rockType == NONE){
        rockNorth = 0.;
    }else if(m_rockType == UPDOWN){
        if(m_DECZenith <= 0) rockNorth *= -1.;
    }else if(m_rockType == SLEWING){
        //slewing is experimental
        if(m_DECZenith <= 0) rockNorth *= -1.;
        if(m_DECZenith <= 10 || m_DECZenith <= 10){
            rockNorth -= rockNorth*(m_DECZenith/10.);
        }
    }else if(m_rockType == ONEPERORBIT){
        //this needs an implementation - it only rocks one way now!
    }
    // now, we want to find the proper transformation for the rocking angles:
    //HepRotation rockRot(dirZ.dir() , -inclination*cos(orbitPhase));
    //HepRotation rockRot2(dirX.dir() , rockNorth);
    HepRotation rockRot;
    rockRot/*.rotateZ(inclination*cos(orbitPhase))*/.rotateX(rockNorth);

    return rockRot;
}

/*Rotation GPS::CELTransform(double seconds){
    // Purpose:  Return the 3x3 matrix which transforms a vector from a local 
    // coordinate system to a galactic coordinate system.
    double time = m_earthOrbit->dateFromSeconds(seconds);
    double degsPerRad = 180./M_PI;
    Rotation gal;//,cel;
    //gal is the matrix which rotates (cartesian)celestial coordiantes into (cartesian)galactic ones
    gal.rotateZ(-282.25/degsPerRad).rotateX(-62.6/degsPerRad).rotateZ(33./degsPerRad);
    
    //cel is the rotation matrix from astro which rotates local coordinates to celestial ones.
    Rotation cel = m_earthOrbit->CelestialToLocal(time).inverse();
    
    //std::cout << "time is " << seconds << std::endl;
    //m_orbit->displayRotation(cel);

    //so gal*cel should be the matrix that makes local coordiates into galactic ones.
    Rotation glstToGal=gal*cel;
    return glstToGal.inverse();

}*/

Rotation GPS::CELTransform(double seconds){
    // Purpose:  Return the 3x3 matrix which transforms a vector from a galactic 
    // coordinate system to a local coordinate system.
    using namespace astro;
    double degsPerRad = 180./M_PI;
    Rotation gal;//,cel;
    double time = m_earthOrbit->dateFromSeconds(seconds);
    
    //double inclination = m_earthOrbit->inclination();
    //double orbitPhase = m_earthOrbit->phase(time);
    m_position = m_earthOrbit->position(time);
    
    //first make the directions for the x and Z axes, as well as the zenith direction.
    //before rotation, the z axis points along the zenith:
    SkyDir dirZ(m_position.unit());
    SkyDir dirX(dirZ.ra()-90., 0.0);
    
    //rotate the x direction so that the x direction points along the orbital direction.
    //dirX().rotate(dirZ.dir() , inclination*cos(orbitPhase));

    //so now we know where the x and z axes of the zenith-pointing frame point in the celestial frame.
    //what we want now is to make cel, where
    //cel is the matrix which rotates (cartesian)local coordinates into (cartesian)celestial ones
    Rotation cel(dirX() , dirZ().cross(dirX()) , dirZ());

    //std::cout << "time is " << seconds << std::endl;
    //m_orbit->displayRotation(cel);

    //gal is the matrix which rotates (cartesian)celestial coordiantes into (cartesian)galactic ones
    gal.rotateZ(-282.25/degsPerRad).rotateX(-62.6/degsPerRad).rotateZ(33./degsPerRad);
    //so gal*cel should be the matrix that makes local coordiates into galactic ones.
    Rotation glstToGal=gal*cel;
    return glstToGal.inverse();

}

Rotation GPS::transformCelToGlast(double seconds){
    // Purpose:  Return the 3x3 matrix which transforms a vector from a celestial 
    // coordinate system (like a SkyDir vector) to a local coordinate system (like the FluxSvc output).
    using namespace astro;
    double degsPerRad = 180./M_PI;

    double time = m_earthOrbit->dateFromSeconds(seconds);
    
    //double inclination = m_earthOrbit->inclination();
    //double orbitPhase = m_earthOrbit->phase(time);
    m_position = m_earthOrbit->position(time);
    
    //first make the directions for the x and Z axes, as well as the zenith direction.
    //before rotation, the z axis points along the zenith:
    SkyDir dirZ(m_position.unit());
    SkyDir dirX(dirZ.ra()-90., 0.0);
    
    //rotate the x direction so that the x direction points along the orbital direction.
    //dirX().rotate(dirZ.dir() , inclination*cos(orbitPhase));

    //so now we know where the x and z axes of the zenith-pointing frame point in the celestial frame.
    //what we want now is to make cel, where
    //cel is the matrix which rotates (cartesian)local coordinates into (cartesian)celestial ones
    Rotation cel(dirX() , dirZ().cross(dirX()) , dirZ());
    return cel.inverse();
}

Rotation GPS::transformGlastToGalactic(double seconds){
    return (/*m_orbit->*/CELTransform(seconds).inverse())*(rockingAngleTransform(seconds).inverse());
}

void GPS::getPointingCharacteristics(double seconds){
    //this is being used by exposureAlg right now, and should be reworked
    //to use the rest of this class.
    using namespace astro;
    
    double time = m_earthOrbit->dateFromSeconds(seconds);
    
    double inclination = m_earthOrbit->inclination();
    double orbitPhase = m_earthOrbit->phase(time);
    m_position = m_earthOrbit->position(time);
    
    //first make the directions for the x and Z axes, as well as the zenith direction.
    SkyDir dirZenith(m_position.unit());
    //before rotation, the z axis points along the zenith:
    SkyDir dirZ(m_position.unit());
    SkyDir dirX(dirZ.ra()-90., 0.0);
    
    //rotate the x direction so that the x direction points along the orbital direction.
    dirX().rotate(dirZ.dir() , inclination*cos(orbitPhase));
    
    //now set the zenith direction before the rocking.
    m_RAZenith = dirZ.ra();
    m_DECZenith = dirZ.dec();
    
    // now, we want to find the proper transformation for the rocking angles:
    //HepRotation rockRot(Hep3Vector(0,0,1).cross(dirZ.dir()) , rockNorth);    
    //and apply the transformation to dirZ and dirX:
    double rockNorth = m_rockDegrees*M_PI/180;
    if(m_DECZenith <= 0) rockNorth *= -1.;
    
    dirZ().rotate(dirX.dir() , rockNorth);
    //dirX().rotate(dirX.dir() , rockNorth*M_PI/180.);//unnecessary
    
    m_RAX = dirX.ra();
    m_RAZ = dirZ.ra();
    m_DECX = dirX.dec();
    m_DECZ = dirZ.dec();
    
    //a test - to ensure the rotation went properly
    //std::cout << " degrees between xhat and zhat directions: " <<
    //    dirZ.difference(dirX)*180./M_PI << std::endl;
}

void GPS::setRockType(int rockType){
    m_rockType = NONE;
    if(rockType == 1) m_rockType = UPDOWN;
    if(rockType == 2) m_rockType = SLEWING;
    if(rockType == 3) m_rockType = ONEPERORBIT;
}
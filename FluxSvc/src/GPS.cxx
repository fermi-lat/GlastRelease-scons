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
    m_orbit->setpitch(p);
}

void GPS::yaw ( double y )
{
    m_orbit->setyaw(y);
}

void GPS::roll ( double r )
{
    m_orbit->setroll(r);
}

void GPS::lat ( double l )
{
    m_orbit->setlatitude(l);
}

void GPS::lon ( double l )
{
    m_orbit->setlongitude(l);
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
    m_orbit->setlatitude(c.lat());
    m_orbit->setlongitude(c.lon());
    m_orbit->setpitch(c.pitch());
    m_orbit->setyaw(c.yaw());
    m_orbit->setroll(c.roll());
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


Vector GPS::earthToGlast(Vector launchDir){
    Vector retlaunch=launchDir;
    double theta=rotateAngles().first;
    double phi=rotateAngles().second;
    retlaunch.rotateX(theta).rotateZ(phi);
    return retlaunch;
}

Vector GPS::galaxyToGlast(Vector launchDir){
    
    //code for changing vector into theta, phi (actually l,b, as these are the declared coordinates for the source)
    double lpri = (-launchDir).theta();
    double bpri = (-launchDir).phi();
    
    //then get galactic coordinates of source corresponding to these...
    //double lpri=thetapri*cos(phipri);
    //double bpri=thetapri*sin(phipri);
    
    //now we want the (l.g) coordinates of the GLAST satellite.
    std::pair<double,double> glastpoint = galPositionOfGlast();
    double lg=glastpoint.first;
    double bg=glastpoint.second;
    
    //then get glast-relative theta and phi..
    double thetat=sqrt(pow(lg-lpri,2)+pow(bg-bpri,2));
    double phit=0.; ///FIX THIS!!!
    
    //then turn the relative theta, phi back into a vector...
    
    Vector transformDir(sin(thetat)*cos(phit),sin(thetat)*sin(phit),cos(thetat));
    return transformDir;
}



std::pair<double,double> GPS::galPositionOfGlast(){
    
    //lat, lon, sidereal??? time becomes J2000.0 coordinates here:
    double ra=((orbittime()/(24.0*60.0*60.0))*360.0)+lon();
    if(ra > 360.0) ra-=360.0;
    if(ra < 0.0 ) ra+=360.0;
    double dec=lat();
    
    //can consider GLAST to be zenith-pointing here....m_rotangles takes care of the rest...
    double pi = 3.141592653589793238;
    //double ra=0;
    //double dec=0;
    double b,l,otherl,cosl,cosb,sinb,sinl,cosra,sinra,cosdec,sindec;
    double radperdeg=(2.0*pi)/360.0;
    cosra=cos(ra*radperdeg);
    sinra=sin(ra*radperdeg);
    
    cosdec=cos(dec*radperdeg);
    sindec=sin(dec*radperdeg);
    double cosbcosl=(cosdec*cosra*-0.0548755)+(cosdec*sinra*-0.87343711)+(sindec*-0.4838349858);
    double cosbsinl=(cosdec*cosra*0.49410945)+(cosdec*sinra*-0.444829589)+(sindec*0.7469822518);
    sinb=(cosdec*cosra*-0.8676661358)+(cosdec*sinra*-0.198076386)+(sindec*0.4559837957);
    b=asin(sinb); //BUT WHAT OF b < -90 or > 90?
    cosb=cos(b);
    cosl=cosbcosl/cosb;
    sinl=cosbsinl/cosb;
    otherl=acos(cosl);
    l=asin(sinl);
    //if(l - otherl >= 0.001) std::cout << "unequal stuff! " << std::endl;
    b=b/radperdeg;
    l=l/radperdeg;

    //return the galactic-pointing coordinates of glast: (b,l)
    return std::make_pair<double,double>(b,l);
}





//access m_rotangles
std::pair<double,double> GPS::rotateAngles(){
    return m_rotangles;
    
}

//set m_rotangles
void GPS::rotateAngles(std::pair<double,double> coords){
    m_rotangles=coords;
}


Rotation GPS::rockingAngleTransform(double time){
    
    Rotation gal;   
    //and here we construct the rotation matrix
    double zenithPhase = m_rotangles.first;
    double offZenith = m_rotangles.second;
    //gal.rotateZ(zenithPhase).rotateX(offZenith);
    gal.rotateX(offZenith).rotateZ(zenithPhase);

    return gal;
}

Rotation GPS::transformGlastToGalactic(double time){
    return (m_orbit->CELtransform(time).inverse())*(rockingAngleTransform(time).inverse());
}


// $Header$
#include "Event/MonteCarlo/D2Entry.h"

using namespace Event;

void D2Entry::writeOut(MsgStream& log) const
{
    log << MSG::INFO << "livetime =" << livetime() << endreq;
    log << MSG::INFO << "rax =" << RAX() << endreq;
    log << MSG::INFO << "decx =" << DECX() << endreq;
    log << MSG::INFO << "raz =" << RAZ() << endreq;
    log << MSG::INFO << "decz =" << DECZ() << endreq;
    log << MSG::INFO << "lat =" << lat() << endreq;
    log << MSG::INFO << "lon =" << lon() << endreq;
    log << MSG::INFO << "alt =" << alt() << endreq;
    log << MSG::INFO << "ramoon =" << RAMoon() << endreq;
    log << MSG::INFO << "decmoon =" << DECMoon() << endreq;
    log << MSG::INFO << "rasun =" << RASun() << endreq;
    log << MSG::INFO << "decsun =" << DECSun() << endreq;
    log << MSG::INFO << "Interval Start = " << intrvalstart() << endreq;
    log << MSG::INFO << "Interval End = " << intrvalend() << endreq;
    log << MSG::INFO << "Position[x] = " << posX() << endreq;
    log << MSG::INFO << "Position[y] = " << posY() << endreq;
    log << MSG::INFO << "Position[z] = " << posZ() << endreq;
    log << MSG::INFO << "RAZenith = " << RAZenith() << endreq;
    log << MSG::INFO << "DECZenith = " << DECZenith() << endreq;
    log << MSG::INFO << "In SAA ? = " << SAA() << endreq;
    return;
}

//! completely initialize a new object.
void D2Entry::init(double  posX,double  posY,double  posZ,double RAX,double RAZ,double DECX,double DECZ,double  RAZenith,double  DECZenith, double lat,double lon,double alt,
                       double  intrvalstart,double  intrvalend,double livetime, double RAMoon, double DECMoon,double RASun,double DECSun,bool  SAA){
    
    //set all the object parameters.
    m_RAX = RAX;
    m_RAZ = RAZ;
    m_DECX = DECX;
    m_DECZ = DECZ;
    m_lat = lat;
    m_lon = lon;
    m_alt = alt;
    m_livetime = livetime;
    m_RAMoon = RAMoon;
    m_DECMoon = DECMoon;
    m_RASun = RASun;
    m_DECSun = DECSun;
    m_intrvalstart = intrvalstart;
    m_intrvalend = intrvalend;
    m_posX = posX;
    m_posY = posY;
    m_posZ = posZ;
    m_RAZenith = RAZenith;
    m_DECZenith = DECZenith;
    m_SAA = SAA;
    
}

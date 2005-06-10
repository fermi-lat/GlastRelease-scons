#include "Event/Recon/CalRecon/CalMipClasses.h"

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
void Event::CalMipXtal::initialize(Event::CalXtalRecData* xtalData, double d2C, bool free)
{
    m_xtalData = xtalData;
    m_d2C      = d2C;
    m_free     = free;
}

//-----------------------------------------------------------------------------------------------------------------
 
void Event::CalMipXtal::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive() ) 
    {
        log << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endreq;
        log << " Print in CalMipXtal s Class"             << endreq;
     
        log << "Free=" << m_free                          << endreq;
        log << "D2C =" << m_d2C                           << endreq;
        log << "Ener=" << m_xtalData->getEnergy()         << endreq;
        log << "PosX=" << m_xtalData->getPosition().x()   << endreq;
        log << "PosY=" << m_xtalData->getPosition().y()   << endreq;
        log << "PosZ=" << m_xtalData->getPosition().z()   << endreq;
     
        log << " Print End"                               << endreq;
    }
    log << endreq;
}

std::ostream& Event::CalMipXtal::fillStream( std::ostream& s ) const 
{ 
    s << "Free=" << m_free                          << "\n"
      << "D2C =" << m_d2C                           << "\n"
      << "Ener=" << m_xtalData->getEnergy()         << "\n"
      << "PosX=" << m_xtalData->getPosition().x()   << "\n"
      << "PosY=" << m_xtalData->getPosition().y()   << "\n"
      << "PosZ=" << m_xtalData->getPosition().z();
  
  return s; 
}

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
void Event::CalMipTrack::initialize(Point point, Vector vector, int ndofTrack, double ki2Track, CalMipXtalVec calMipXtalTrack)
{
    m_point     = point;
    m_vector    = vector;
    m_ndofTrack = ndofTrack;
    m_ki2Track  = ki2Track;

    (CalMipXtalVec)(*this) = calMipXtalTrack;
}
//-----------------------------------------------------------------------------------------------------------------
void Event::CalMipTrack::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive() ) 
    {
        log << "-----------------------------------------------------------------" << endreq;
        log << " Print in CalMipTrack s Class" << endreq;

        log <<  "Nh  =" << size()          << endreq;
        log <<  "Ndof=" << m_ndofTrack     << endreq;
        log <<  "Ki2 =" << m_ki2Track      << endreq;
        log <<  "PosX=" << getPoint().x()  << endreq;
        log <<  "PosY=" << getPoint().y()  << endreq;
        log <<  "PosZ=" << getPoint().z()  << endreq;
        log <<  "DirX=" << getDir().x()    << endreq;
        log <<  "DirY=" << getDir().y()    << endreq;
        log <<  "DirZ=" << getDir().z()    << endreq;

        int comptor=0;

        for(Event::CalMipXtalVec::const_iterator xTalIter = begin(); xTalIter != end(); xTalIter++)
        {
            log << "-------------------------------------------" << endreq;
            log << "Xtal No="<<comptor                           << endreq;
            log << "---------------"                             << endreq;
            comptor++;
            Event::CalMipXtal calMipXtal=*xTalIter;
            calMipXtal.writeOut(log);
        }
 
        log << " Print - End"  << endreq;
    }
    log << endreq;
}
//-----------------------------------------------------------------------------------------------------------------
std::ostream& Event::CalMipTrack::fillStream( std::ostream& s ) const
{
    s << "-----------------------------------------------------------------" << "\n"
      << " Print in CalMipTrack s Class" << "\n"

      <<  "Nh  =" << size()          << "\n"
      <<  "Ndof=" << m_ndofTrack     << "\n"
      <<  "Ki2 =" << m_ki2Track      << "\n"
      <<  "PosX=" << getPoint().x()  << "\n"
      <<  "PosY=" << getPoint().y()  << "\n"
      <<  "PosZ=" << getPoint().z()  << "\n"
      <<  "DirX=" << getDir().x()    << "\n"
      <<  "DirY=" << getDir().y()    << "\n"
      <<  "DirZ=" << getDir().z();

    int comptor=0;

    for(Event::CalMipXtalVec::const_iterator xTalIter = begin(); xTalIter != end(); xTalIter++)
    {
        s << "-------------------------------------------" << "\n"
          << "Xtal No="<<comptor                           << "\n"
          << "---------------";
        comptor++;
        Event::CalMipXtal calMipXtal=*xTalIter;
        calMipXtal.fillStream(s);
    }

    return s;
}

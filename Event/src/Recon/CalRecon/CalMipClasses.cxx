#include "Event/Recon/CalRecon/CalMipClasses.h"

//-----------------------------------------------------------------------------------------------------------------
void Event::CalMipXtal::initialize(Event::CalXtalRecData* xtalData, double d2C, bool free, bool freeC0, double ecor)
{
    m_xtalData = xtalData;
    m_d2C      = d2C;
    m_free     = free;
    m_freeC0   = freeC0;
    m_ecor     = ecor;
}

//-----------------------------------------------------------------------------------------------------------------
void Event::CalMipXtal::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive() ) 
    {
        int it=m_xtalData->getPackedId().getTower();
        int il=m_xtalData->getPackedId().getLayer();
        int ic=m_xtalData->getPackedId().getColumn();
        log << "---> writeOut CalMipXtal tow / lay / col = " << it << " " << il << " " << ic << endreq;
        log << "---> x / y /z = " << m_xtalData->getPosition().x() << " " << m_xtalData->getPosition().y() << " " << m_xtalData->getPosition().z() << endreq;
        log << "---> free / freeC0 / d2C / e / ecor = " << m_free << " " << m_freeC0  << " " << m_d2C << " " << m_xtalData->getEnergy() << " " << m_ecor << endreq;
    }
}

//-----------------------------------------------------------------------------------------------------------------
std::ostream& Event::CalMipXtal::fillStream( std::ostream& s ) const 
{ 
    s << "D2C ="   << m_d2C      << "\n"
      << "Free="   << m_free     << "\n"
      << "FreeC0=" << m_freeC0   << "\n"
      << "eocr="   << m_ecor     << "\n";
  return s; 
}

//-----------------------------------------------------------------------------------------------------------------
void Event::CalMipTrack::initialize(CalMipXtalVec calMipXtalTrack, Point point, Vector dir, double dt2C, double d2Edge, int calEdge, double arcLen, double ecor, double ecorRms, double chi2, double erm)
{
  m_point    = point;
  m_dir      = dir;
  m_dt2C     = dt2C;
  m_d2Edge   = d2Edge;
  m_calEdge  = calEdge;
  m_arcLen   = arcLen;
  m_ecor     = ecor;
  m_ecorRms  = ecorRms;
  m_chi2     = chi2;
  m_erm      = erm;

  (CalMipXtalVec)(*this) = calMipXtalTrack;
}

//-----------------------------------------------------------------------------------------------------------------
void Event::CalMipTrack::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive() ) 
    {
      log << "------------------------------------------------------------" << endreq;
      log << "---> writeOutPrint calMipTrack size = " << size() << endreq;
      log << "---> point =" << getPoint().x() << " " << getPoint().y() << " " << getPoint().z() << endreq;
      log << "---> dir   =" << getDir().x()   << " " << getDir().y()   << " " << getDir().z()   << endreq;
      log << "---> dt2C / d2Edge / calEdge / arcLen = " << m_dt2C << " " << m_d2Edge << " " << m_calEdge << " "<< m_arcLen << endreq;
      log << "---> ecor / ecorRms / chi2 / erm = " << m_ecor << " " << m_ecorRms << " " << m_chi2 <<" " << m_erm << endreq;
      int counter=0;
      for(Event::CalMipXtalVec::const_iterator xTalIter = begin(); xTalIter != end(); xTalIter++)
      {
        log << "---> Xtal No="<< counter << endreq;
        counter++;
        Event::CalMipXtal calMipXtal=*xTalIter;
        calMipXtal.writeOut(log);
      }
      log << "------------------------------------------------------------" << endreq;
    }
}

//-----------------------------------------------------------------------------------------------------------------
std::ostream& Event::CalMipTrack::fillStream( std::ostream& s ) const
{
    s << "---> writeOutPrint calMipTrack size = " << getNh() << "\n"
      << "---> point =" << getPoint().x() << " " << getPoint().y() << " " << getPoint().z() << "\n"
      << "---> dir   =" << getDir().x()   << " " << getDir().y()   << " " << getDir().z()   << "\n"
      << "---> dt2C / d2Edge / calEdge / arcLen = " << m_dt2C << " " << m_d2Edge << " " << m_calEdge << " "<< m_arcLen << "\n"
      << "---> ecor / ecorRms / chi2 / erm = " << m_ecor << " " << m_ecorRms << " " << m_chi2 << " " << m_erm;

    int counter=0;

    for(Event::CalMipXtalVec::const_iterator xTalIter = begin(); xTalIter != end(); xTalIter++)
    {
        s << "---> Xtal No="<< counter << endreq;
        counter++;
        Event::CalMipXtal calMipXtal=*xTalIter;
        calMipXtal.fillStream(s);
    }

    return s;
}

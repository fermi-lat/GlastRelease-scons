
#include "Event/Recon/TkrRecon/TkrVertex.h"

using namespace Event;

TkrVertex::TkrVertex(idents::TkrId tkrID, double energy, double quality, double chisq, 
			       double radlen, double doca, double s1, double s2, double z,
				   TkrTrackParams params):
           m_params(params), m_energy(energy), m_quality(quality), m_chiSquare(chisq), m_statusBits(0), 
		   m_doca(doca), m_arcLen1(s1), m_arcLen2(s2), m_radlen(radlen), m_vtxID(tkrID)
{
    m_position   = Point(params(1), params(3), z);
    m_direction  = 	Vector(-params(2), -params(4), -1.).unit();
    m_tracks.clear();
}


void TkrVertex::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive() ) {
        log << " --- TkrVertex::writeOut --- " << endreq
            << " Position      = " << getPosition().x()  << " " 
            << getPosition().y()  << " " << getPosition().z() << endreq
            << " Direction     = " << getDirection().x() << " " 
            << getDirection().y() << " " << getDirection().z() << endreq
            << " Energy        = " << getEnergy();
    }
    log << endreq;
}

std::ostream& TkrVertex::fillStream( std::ostream& s ) const 
{ 
  s << " Position      = " << getPosition().x()  << " " << getPosition().y()  << " " << getPosition().z()  << "\n"
    << " Direction     = " << getDirection().x() << " " << getDirection().y() << " " << getDirection().z() << "\n"
    << " Energy        = " << getEnergy() << "\n"
    << " quality       = " << getQuality();
  
  return s; 
}

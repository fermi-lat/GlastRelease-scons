
#include "Event/Recon/TkrRecon/TkrVertex.h"

using namespace Event;

TkrVertex::TkrVertex(int ilyr, int itwr, double energy, double quality, const Ray& testRay)
{
    m_position   = testRay.position();
    m_direction  = testRay.direction();
    m_energy     = energy;
    m_quality    = quality;
    m_firstLayer = ilyr;
    m_itower     = itwr;

    if (m_direction.mag() != 1.) 
    {
        m_direction = m_direction.unit();
    }

    m_vertexPar  = TkrFitPar(m_position.x(),m_direction.x()/m_direction.z(),m_position.y(),m_direction.y()/m_direction.z());
    m_vertexCov  = TkrFitMatrix();

    m_vertexCov(1,1) = 1.;
    m_vertexCov(2,2) = 1.;
    m_vertexCov(3,3) = 1.;
    m_vertexCov(4,4) = 1.;

    m_tracks.clear();
}


void TkrVertex::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrVertex::writeOut --- " << endreq;

    log << MSG::DEBUG << " Position      = " << getPosition().x()  << " " << getPosition().y()  << " " << getPosition().z() << endreq;
    log << MSG::DEBUG << " Direction     = " << getDirection().x() << " " << getDirection().y() << " " << getDirection().z() << endreq;
    log << MSG::DEBUG << " Energy        = " << getEnergy() << endreq;
    log << MSG::DEBUG << " first Layer   = " << getLayer() << endreq;
    log << MSG::DEBUG << " Tower         = " << getTower() << endreq;
}

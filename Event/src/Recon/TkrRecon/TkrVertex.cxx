
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
    log << MSG::DEBUG;
    if (log.isActive() ) {
        log << " --- TkrVertex::writeOut --- " << endreq
            << " Position      = " << getPosition().x()  << " " 
            << getPosition().y()  << " " << getPosition().z() << endreq
            << " Direction     = " << getDirection().x() << " " 
            << getDirection().y() << " " << getDirection().z() << endreq
            << " Energy        = " << getEnergy() << endreq
            << " first Layer   = " << getLayer() << endreq
            << " Tower         = " << getTower();
    }
    log << endreq;
}

double        TkrVertex::getQuality()            const 
{
    return m_quality;
};
double        TkrVertex::getEnergy(TrackEnd)     const 
{
    return m_energy;
}
int           TkrVertex::getLayer(TrackEnd )     const 
{
    return m_firstLayer;
}
int           TkrVertex::getTower(TrackEnd )     const 
{
    return m_itower;
}
Point         TkrVertex::getPosition(TrackEnd )  const 
{
    return m_position;
}
Vector        TkrVertex::getDirection(TrackEnd ) const 
{
    return m_direction;
}
Ray           TkrVertex::getRay(TrackEnd )       const 
{
    return Ray(getPosition(),getDirection());
}
TkrFitPar     TkrVertex::getTrackPar(TrackEnd )  const 
{
    return m_vertexPar;
}
double        TkrVertex::getTrackParZ(TrackEnd ) const 
{
    return m_position.z();
}
TkrFitMatrix  TkrVertex::getTrackCov(TrackEnd )  const 
{
    return m_vertexCov;}
bool          TkrVertex::empty(int)              const 
{
    return m_firstLayer >= 0;
}


/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include <algorithm>
#include "Event/Recon/TkrRecon/TkrPatCand.h"

using namespace Event;

TkrPatCand::TkrPatCand(int layer, int tower, double energy, double quality, const Ray& testRay)
{
    //Zero out the candidate hit vector
    m_hits.clear();

    //Initialize the TkrBase info
    m_position   = testRay.position();
    m_direction  = testRay.direction();
    m_firstLayer = layer;
    m_itower     = tower;
    m_energy     = energy;
    m_quality    = quality;

    return;
}

TkrPatCand::~TkrPatCand()
{
}

void TkrPatCand::addCandHit(TkrCluster* pCluster)
{
    m_hits.push_back(TkrPatCandHit(pCluster));
    std::sort(m_hits.begin(),m_hits.end());

    return;
}

void TkrPatCand::addCandHit(TkrPatCandHit CandHit)
{
    m_hits.push_back(CandHit);
    std::sort(m_hits.begin(),m_hits.end());

    return;
}

int TkrPatCand::lastLayer()
{
    int finalLayer   = getLayer();

    if (m_hits.size() > 0)
    {
        finalLayer = getFoLPlane(TkrRecInfo::End)->PlaneIndex();
    }

    return finalLayer;
}

TkrPatCandHit* TkrPatCand::getFoLPlane(TrackEnd end)
{
    if (m_hits.size() == 0) 
    {
        return 0;
    }
    else
    {
        if (end == TkrRecInfo::Start) return &m_hits.front();
        else                          return &m_hits.back();
    }
}

TkrFitPar TkrPatCand::getTrackPar(TrackEnd end)    const
{
    return TkrFitPar(m_position.x(),m_direction.x(),m_position.y(),m_direction.y());
}
    
TkrFitMatrix TkrPatCand::getTrackCov(TrackEnd end)    const
{
    TkrFitMatrix cov;

    cov(1,1) = 1.;
    cov(2,2) = 1.;
    cov(3,3) = 1.;
    cov(4,4) = 1.;

    return cov;
}

void TkrPatCand::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrPatCandHit::writeOut --- " << endreq;
    log << MSG::DEBUG << " Position      = " << getPosition().x()  << " " << getPosition().y()  << " " << getPosition().z()  << endreq;
    log << MSG::DEBUG << " Direction     = " << getDirection().x() << " " << getDirection().y() << " " << getDirection().z() << endreq;
    log << MSG::DEBUG << " Energy        = " << getEnergy() << endreq;
    log << MSG::DEBUG << " first Layer   = " << getLayer() << endreq;
    log << MSG::DEBUG << " Tower         = " << getTower() << endreq;
}


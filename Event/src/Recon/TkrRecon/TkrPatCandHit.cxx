/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include "Event/Recon/TkrRecon/TkrPatCandHit.h"

using namespace Event;

TkrPatCandHit::TkrPatCandHit(TkrCluster* pCluster)
{
    m_position   = pCluster->position();
    m_hitIndex   = pCluster->id();
    m_towerIndex = pCluster->tower();
    m_planeIndex = pCluster->getTkrId().getPlane();
    m_view       = pCluster->getTkrId().getView();

    return;
}

TkrPatCandHit::TkrPatCandHit(unsigned int hitId, const Point& pos, unsigned int tower, 
                             unsigned int layer, int v) {

    m_position = pos;
    m_hitIndex = hitId;
    m_towerIndex = tower;
    m_planeIndex = layer;
    m_view = v;
    return;
}


void TkrPatCandHit::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive()) {
        log << " --- TkrPatCandHit::writeOut --- " << endreq
            << " Position= " << Position().x() << ", " 
            <<Position().y() << ", " << Position().z() << endreq
            << " Tower: " << TowerIndex() << ", Layer: " 
            << PlaneIndex() << ", view: " << View();
    }
    log << endreq;
}


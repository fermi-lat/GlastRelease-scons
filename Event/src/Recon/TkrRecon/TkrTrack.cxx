//------------------------------------------------------------------------------
//
//     TkrTrack
//
//     Implementation of the Tracker Fit Track output class
//				
//------------------------------------------------------------------------------


#include "Event/Recon/TkrRecon/TkrTrack.h"

using namespace Event;

// Constructor takes no arguments and basically just initializes to
// a null state. In the class' current incarnation, it is expected to
// be inherited from a class which actually does the track fit. 
TkrTrack::TkrTrack()
{
    m_statusBits       = 0;

    m_initialEnergy    = 0.;
    m_initialPosition  = Point(0.,0.,0.);
    m_initialDirection = Vector(0.,0.,0.);

    m_chiSquareFilter  = 0.;
    m_chiSquareSmooth  = 0.;
    m_rmsResid         = 0.;
    m_Quality          = 0.;

    m_KalmanEnergy     = 0.;
    m_KalmanThetaMS    = 0.;
    m_Xgaps            = 0;
    m_Ygaps            = 0;
    m_XistGaps         = 0;
    m_YistGaps         = 0;

    m_numSegmentPoints = 0;
    m_chisqSegment     = 0.;
    m_nxHits           = 0;
    m_nyHits           = 0;
    m_KalmanEnergyErr  = 0.;
    m_TkrCal_radlen    = 0.; 


    clear();
}

TkrTrack::~TkrTrack()
{
    // Is this track registered in the TDS?
    if (parent() == 0)
    {
        // If not then probably the TkrTrackHits aren't either in which case we need
        // to take care of deleting them
        for(SmartRefVector<Event::TkrTrackHit>::iterator hitIter  = begin();
                                                         hitIter != end(); 
                                                         hitIter++)
        {
            Event::TkrTrackHit* trackHit = *hitIter;

            if (trackHit->parent() == 0) delete trackHit;
        }
    }

    clear();
}

void TkrTrack::writeOut(MsgStream& log) const
{
    
    log << MSG::DEBUG;
    if (log.isActive() ) {
//        log << " --- TkrTrack::writeOut --- "            << endreq
//            << " Position      = " << getPosition().x()  << " " 
//            << getPosition().y()  << " " << getPosition().z()  << endreq
//            << " Direction     = " << getDirection().x() << " " 
//            << getDirection().y() << " " << getDirection().z() << endreq
//            << " Energy        = " << getEnergy() << endreq
//            << " first Layer   = " << getLayer() << endreq
//            << " Tower         = " << getTower() << endreq
//            << " quality       = " << getQuality()       << endreq
//            << " num m_hits    = " << getNumHits();
    }
    log << endreq;
}


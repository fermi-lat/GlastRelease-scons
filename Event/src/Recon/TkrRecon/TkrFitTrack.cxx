//------------------------------------------------------------------------------
//
//     TkrFitTrack
//
//     Implementation of the Tracker Fit Track output class
//				
//------------------------------------------------------------------------------


#include "Event/Recon/TkrRecon/TkrFitTrack.h"

using namespace Event;

// Constructor takes no arguments and basically just initializes to
// a null state. In the class' current incarnation, it is expected to
// be inherited from a class which actually does the track fit. 
TkrFitTrack::TkrFitTrack()
{ 
    m_chisq       = -1e6;
    m_chisqSmooth = -1e6;
    m_rmsResid    = 0.;
    m_Q           = -1e6;
    m_KalEnergy   = 0.;
    m_KalThetaMS  = 0.;
    m_Xgaps       = 0;
    m_Ygaps       = 0;
    m_XistGaps    = 0;
    m_YistGaps    = 0;

    m_hits.clear();
}


void TkrFitTrack::initializeInfo(unsigned int xgaps, unsigned int ygaps, 
                                 unsigned int x1st, unsigned int y1st) {
    m_Xgaps    = xgaps;
    m_Ygaps    = ygaps;
    m_XistGaps = x1st;
    m_YistGaps = y1st;

}

void TkrFitTrack::initializeQual(double chiSq, double ChiSqSmooth, double rms, double quality, double e, double ms)
{
    m_chisq       = chiSq;
    m_chisqSmooth = ChiSqSmooth;
    m_rmsResid    = rms;
    m_Q           = quality;
    m_KalEnergy   = e;
    m_KalThetaMS  = ms;
}


TkrFitTrack::~TkrFitTrack()
{
    clear();
}

bool TkrFitTrack::empty(int numHits) const
{
    bool empty = false;
    if (getLayer()     < 0)       empty = true;
    if (getNumHits()   < numHits) empty = true;
    if (getChiSquare() < 0.)      empty = true;

    return empty;
}

void TkrFitTrack::clear()
{   
    m_hits.clear();
    
    m_Q            = -1e6;
}

Vector TkrFitTrack::getDirection(TrackEnd end) const
{
    TkrFitPar trkPar = getFoLPlane(end).getHit(TkrFitHit::SMOOTH).getPar();

    return Vector(-trkPar.getXSlope(),-trkPar.getYSlope(),-1.).unit();
}



void TkrFitTrack::writeOut(MsgStream& log) const
{
    
    log << MSG::DEBUG << " --- TkrFitTrack::writeOut --- "            << endreq;
    log << MSG::DEBUG << " Position      = " << getPosition().x()  << " " << getPosition().y()  << " " << getPosition().z()  << endreq;
    log << MSG::DEBUG << " Direction     = " << getDirection().x() << " " << getDirection().y() << " " << getDirection().z() << endreq;
    log << MSG::DEBUG << " Energy        = " << getEnergy() << endreq;
    log << MSG::DEBUG << " first Layer   = " << getLayer() << endreq;
    log << MSG::DEBUG << " Tower         = " << getTower() << endreq;
    log << MSG::DEBUG << " quality       = " << getQuality()       << endreq;
    log << MSG::DEBUG << " num m_hits    = " << getNumHits()       << endreq;
}


TkrFitPlane TkrFitTrack::getFoLPlane(TrackEnd end) const
{
    if (m_hits.size() == 0) 
    {
        return TkrFitPlane();
    }
    else
    {
        if (end == TkrFitTrackBase::Start) return m_hits.front();
        else                               return m_hits.back();
    }
}

Ray TkrFitTrack::getRay(TrackEnd) const 
{
    return Ray(getPosition(),getDirection());
}

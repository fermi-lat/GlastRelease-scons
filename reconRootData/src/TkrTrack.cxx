#include "reconRootData/TkrTrack.h"

ClassImp(TkrTrack)

// Default constructor
TkrTrack::TkrTrack() 
{
    // Make sure everything is properly zeroed
    m_id          = -9999;
    m_ChiSq       = -9999.;
    m_ChiSqSmooth = -9999.;
    m_rmsResid    = 0.;
    m_Q           = 0.;
    m_KalEnergy   = 0.;
    m_KalThetaMS  = 0.;
    m_Xgaps       = 0;
    m_Ygaps       = 0;
    m_XistGaps    = 0;
    m_YistGaps    = 0;

    m_hits.clear();
}

void TkrTrack::initializeInfo(UInt_t id, UInt_t xgaps, UInt_t ygaps, UInt_t x1st, UInt_t y1st)
{
    m_id       = id;
    m_Xgaps    = xgaps;
    m_Ygaps    = ygaps;
    m_XistGaps = x1st;
    m_YistGaps = y1st;

    m_hits.clear();
}

void TkrTrack::initializeQaul(Double_t chiSq, Double_t ChiSqSmooth, Double_t rms, Double_t quality, Double_t e, Double_t ms)
{
    m_ChiSq       = chiSq;
    m_ChiSqSmooth = ChiSqSmooth;
    m_rmsResid    = rms;
    m_Q           = quality;
    m_KalEnergy   = e;
    m_KalThetaMS  = ms;
}

const TkrHitPlane& TkrTrack::getEndHit(TKREND end)
{
    // If no hits put a bogus one in the list so we don't blow up
    if (m_hits.size() == 0) m_hits.push_back(TkrHitPlane());

    // Return the TkrHitPlane object at the desired end of the track
    if (end == TKREND::FIRSTHIT) return m_hits.front();
    else                         return m_hits.back();
}

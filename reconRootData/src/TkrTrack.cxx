#include "reconRootData/TkrTrack.h"

ClassImp(TkrTrack)

TkrTrack::TkrTrack() {
    Clear();
}

TkrTrack::~TkrTrack() {
    Clear();
}

void TkrTrack::Clear(Option_t *option) {
    // Make sure everything is properly zeroed
    m_id          = -9999;
    m_chiSq       = -9999.;
    m_chiSqSmooth = -9999.;
    m_rmsResid    = 0.;
    m_q           = 0.;
    m_kalEnergy   = 0.;
    m_kalThetaMS  = 0.;
    m_xGaps       = 0;
    m_yGaps       = 0;
    m_XistGaps    = 0;
    m_YistGaps    = 0;

    m_hits.clear();
}

void TkrTrack::Print(Option_t *option) const {
    TObject::Print(option);
}


void TkrTrack::initializeInfo(UInt_t id, UInt_t xgaps, UInt_t ygaps, UInt_t x1st, UInt_t y1st)
{
    m_id       = id;
    m_xGaps    = xgaps;
    m_yGaps    = ygaps;
    m_XistGaps = x1st;
    m_YistGaps = y1st;

    m_hits.clear();
}

void TkrTrack::initializeQual(Double_t chiSq, Double_t ChiSqSmooth, Double_t rms, Double_t quality, Double_t e, Double_t ms)
{
    m_chiSq       = chiSq;
    m_chiSqSmooth = ChiSqSmooth;
    m_rmsResid    = rms;
    m_q           = quality;
    m_kalEnergy   = e;
    m_kalThetaMS  = ms;
}



Double_t TkrTrack::getEnergy(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
       return hitPlane->getEnePlane();
    }
    return 0.0;
}

UInt_t TkrTrack::getLayer(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return hitPlane->getIdPlane();
    }
    return 0;
}

UInt_t TkrTrack::getTower(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return hitPlane->getIdTower();
    }
    return 0;
}

const TkrParams* TkrTrack::getTrackPar(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return &(hitPlane->getHitSmooth().getTkrParams());
    }
    return 0;
}

Double_t TkrTrack::getTrackParZ(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return hitPlane->getZplane();
    } 
    return 0.0;
}

const TkrCovMat* TkrTrack::getTrackCov(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return &(hitPlane->getHitSmooth().getTkrCovMat());
    } 
    return 0;
}

const TkrHitPlane* TkrTrack::getEndHit(TKREND end) const
{
    // If no hits put a bogus one in the list so we don't blow up
    //if (m_hits.size() == 0) m_hits.push_back(TkrHitPlane());
    if (m_hits.size() == 0) return 0;

    // Return the TkrHitPlane object at the desired end of the track
    if (end == TKREND::FIRSTHIT) return &(m_hits.front());
    else                         return &(m_hits.back());
}

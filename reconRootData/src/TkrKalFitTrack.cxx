#include "reconRootData/TkrKalFitTrack.h"

ClassImp(TkrKalFitTrack)

TkrKalFitTrack::TkrKalFitTrack() {
//    m_hits = new TObjArray();
//    m_hits->Clear();
    Clear();
}

TkrKalFitTrack::~TkrKalFitTrack() {
    Clear();
//    delete m_hits;
}

void TkrKalFitTrack::Clear(Option_t *option) {
    // Make sure everything is properly zeroed
    m_status           =  EMPTY;
    m_type             =  0;
    m_id               = -9999;
    m_energy0          = -9999.;
    m_x0               =  TVector3(0.,0.,0.);
    m_dir              =  TVector3(0.,0.,0.);
    m_chisq            = -9999.;
    m_chisqSmooth      = -9999.;
    m_rmsResid         =  0.;
    m_Q                =  0.;
    m_KalEnergy        =  0.;
    m_KalThetaMS       =  0.;
    m_Xgaps            =  0;
    m_Ygaps            =  0;
    m_XistGaps         =  0;
    m_YistGaps         =  0;
    m_numSegmentPoints =  0;
    m_chisqSegment     =  0.;
    m_nxHits           =  0;
    m_nyHits           =  0;
    m_KalEnergyErr     =  0.;
    m_TkrCal_radlen    =  0.;

//    for(int idx=0; idx < m_hits->GetEntries(); idx++) {delete (TkrHitPlane*)m_hits->At(idx);}

    m_hits.clear();
}

void TkrKalFitTrack::Print(Option_t *option) const {
    TObject::Print(option);
}

void TkrKalFitTrack::initializeBase(Status stat, UInt_t type, UInt_t id, Double_t energy0, TVector3& x0, TVector3& dir)
{
    m_status   = stat;
    m_type     = type;
    m_id       = id;
    m_energy0  = energy0;
    m_x0       = x0;
    m_dir      = dir;

    m_hits.clear();
}

void TkrKalFitTrack::initializeQual(Double_t chiSq, Double_t ChiSqSmooth, Double_t rms, Double_t quality, Double_t e, Double_t e_err, Double_t ms)
{
    m_chisq        = chiSq;
    m_chisqSmooth  = ChiSqSmooth;
    m_rmsResid     = rms;
    m_Q            = quality;
    m_KalEnergy    = e;
    m_KalEnergyErr = e_err;
    m_KalThetaMS   = ms;
}

void TkrKalFitTrack::initializeGaps(UInt_t xgaps, UInt_t ygaps, UInt_t x1st, UInt_t y1st)
{
    m_Xgaps    = xgaps;
    m_Ygaps    = ygaps;
    m_XistGaps = x1st;
    m_YistGaps = y1st;
}

void TkrKalFitTrack::initializeKal(UInt_t nSegPoints, Double_t chisqSeg, UInt_t nxHits, UInt_t nyHits,
                                   Double_t radlen)
{
    m_numSegmentPoints = nSegPoints;
    m_chisqSegment     = chisqSeg;
    m_nxHits           = nxHits;
    m_nyHits           = nyHits;
    m_TkrCal_radlen    = radlen;
}


Double_t TkrKalFitTrack::getEnergy(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
       return hitPlane->getEnePlane();
    }
    return 0.0;
}

UInt_t TkrKalFitTrack::getLayer(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return hitPlane->getIdPlane();
    }
    return 0;
}

UInt_t TkrKalFitTrack::getTower(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return hitPlane->getIdTower();
    }
    return 0;
}

const TkrParams* TkrKalFitTrack::getTrackPar(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return &(hitPlane->getHitSmooth().getTkrParams());
    }
    return 0;
}

Double_t TkrKalFitTrack::getTrackParZ(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return hitPlane->getZplane();
    } 
    return 0.0;
}

const TkrCovMat* TkrKalFitTrack::getTrackCov(TKREND end) const {
    const TkrHitPlane *hitPlane = getEndHit(end);
    if (hitPlane) {
        return &(hitPlane->getHitSmooth().getTkrCovMat());
    } 
    return 0;
}

const TkrHitPlane* TkrKalFitTrack::getEndHit(TKREND end) const
{
    // If no hits put a bogus one in the list so we don't blow up
    //if (m_hits->size() == 0) m_hits->push_back(TkrHitPlane());
//    if (m_hits->GetEntries() == 0) return 0;
    if (m_hits.size() == 0) return 0;

    // Return the TkrHitPlane object at the desired end of the track
    if (end == TkrKalFitTrack::FIRSTHIT) return &m_hits.front();
    else                                 return &m_hits.back();
}

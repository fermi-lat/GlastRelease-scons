#include "reconRootData/TkrHitPlane.h"

ClassImp(TkrHitPlane)

TkrHitPlane::TkrHitPlane() {
    Clear();
}

void TkrHitPlane::initializeInfo(UInt_t hit, UInt_t tower, UInt_t plane, AXIS proj,
                                 AXIS nextProj, Double_t z, Double_t energy, Double_t radLen)
{
    m_IdHit      = hit;
    m_IdTower    = tower;
    m_IdPlane    = plane;
    m_Projection = proj;
    m_ProjPlus   = nextProj;
    m_Zplane     = z;
    m_EnePlane   = energy;
    m_RadLen     = radLen;
}

void TkrHitPlane::initializeHits(const TkrFitHit& meas, const TkrFitHit& pred, const TkrFitHit& fit,
                                 const TkrFitHit& smooth, const TkrCovMat& material)
{
    m_HitMeas   = meas;
    m_HitPred   = pred;
    m_HitFit    = fit;
    m_HitSmooth = smooth;

    m_Qmaterial = material;
}

void TkrHitPlane::Clear(Option_t *option) {
    // Make sure everything is properly zeroed
    m_IdHit      = 0;
    m_IdTower    = 0.;
    m_IdPlane    = 0; 
    m_Projection = TkrHitPlane::NONE;
    m_ProjPlus   = TkrHitPlane::NONE;
    m_Zplane     = 0.;
    m_EnePlane   = 0.;
    m_RadLen     = 0.;
}

void TkrHitPlane::Print(Option_t *option) const {
    TObject::Print(option);
}

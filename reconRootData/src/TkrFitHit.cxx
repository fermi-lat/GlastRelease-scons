#include "reconRootData/TkrFitHit.h"

ClassImp(TkrFitHit)

// Default constructor
TkrFitHit::TkrFitHit() 
{
    m_HitType = UNKNOWN;
    m_Params.Zero();
    m_CovMat.Zero();
}

TkrFitHit::TkrFitHit(const TkrFitHit& hit)
{
    m_HitType = hit.getHitType();
    m_Params  = hit.getTkrParams();
    m_CovMat  = hit.getTkrCovMat();
}

void TkrFitHit::initialize(TYPE type, const TkrParams& params, const TkrCovMat& covMat)
{
    m_HitType = type;
    m_Params  = params;
    m_CovMat  = covMat;
}
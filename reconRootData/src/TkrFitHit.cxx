#include "reconRootData/TkrFitHit.h"

ClassImp(TkrFitHit)

TkrFitHit::TkrFitHit() 
{
    m_HitType = UNKNOWN;
}
    
TkrFitHit::TkrFitHit(TYPE type, const TkrParams& params, const TkrCovMat& covMat)
{
    initialize(type, params, covMat);
}

TkrFitHit::TkrFitHit(const TkrFitHit& hit)
{
    m_HitType = hit.getHitType();
    m_Params  = hit.getTkrParams();
    m_CovMat  = hit.getTkrCovMat();
}

TkrFitHit& TkrFitHit::operator=(const TkrFitHit& hit)
{
    m_HitType = hit.getHitType();
    m_Params  = hit.getTkrParams();
    m_CovMat  = hit.getTkrCovMat();

    return *this;
}

void TkrFitHit::initialize(TYPE type, const TkrParams& params, const TkrCovMat& covMat)
{
    m_HitType = type;
    m_Params  = params;
    m_CovMat  = covMat;
}
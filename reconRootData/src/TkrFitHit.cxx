#include "reconRootData/TkrFitHit.h"

ClassImp(TkrFitHit)

// Default constructor
TkrFitHit::TkrFitHit() 
{
    m_HitType = UNKNOWN;
    m_Params.ResizeTo(1);
    m_Params.Zero();
    m_CovMat.ResizeTo(1,1);
    m_Params.Zero();
}
    
TkrFitHit::TkrFitHit(TYPE type, const TkrParams& params, const TkrCovMat& covMat)
{
    initialize(type, params, covMat);
}

TkrFitHit::TkrFitHit(const TkrFitHit& hit)
{
    m_HitType = hit.getHitType();
    m_Params.ResizeTo(hit.getTkrParams());
    m_Params  = hit.getTkrParams();
    m_CovMat.ResizeTo(hit.getTkrCovMat());
    m_CovMat  = hit.getTkrCovMat();
}

TkrFitHit& TkrFitHit::operator=(const TkrFitHit& hit)
{
    m_HitType = hit.getHitType();
    m_Params.ResizeTo(hit.getTkrParams());
    m_Params  = hit.getTkrParams();
    m_CovMat.ResizeTo(hit.getTkrCovMat());
    m_CovMat  = hit.getTkrCovMat();

    return *this;
}

void TkrFitHit::initialize(TYPE type, const TkrParams& params, const TkrCovMat& covMat)
{
    m_HitType = type;
    m_Params.ResizeTo(params);
    m_Params  = params;
    m_CovMat.ResizeTo(covMat);
    m_CovMat  = covMat;
}
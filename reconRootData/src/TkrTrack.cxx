#include "reconRootData/TkrTrack.h"

ClassImp(TkrTrack)

TkrTrack::TkrTrack() 
{
    Clear();
}

TkrTrack::~TkrTrack() 
{
    Clear();
}

void TkrTrack::Clear(Option_t* /* option */) 
{
    m_statusBits         = 0;

    /// The input energy, and current position and direction
    m_initialEnergy      = 0.;
    m_initialPosition    = TVector3(0.,0.,0.);
    m_initialDirection   = TVector3(0.,0.,0.);

    // Track quality information
    m_chiSquareFilter    = 0.;
    m_chiSquareSmooth    = 0.;
    m_nDegreesOfFreedom  = 0;
    m_rmsResid           = 0.;
    m_Quality            = 0.;

    // Information from the fitter
    m_KalmanEnergy       = 0.;
    m_KalmanThetaMS      = 0.;

    // Hit gap information
    m_Xgaps              = 0;
    m_Ygaps              = 0;
    m_XistGaps           = 0;
    m_YistGaps           = 0;

    /// Kalman Filter Track data
    m_numSegmentPoints   = 0;
    m_chisqSegment       = 0.;
    m_nxHits             = 0;
    m_nyHits             = 0;
    m_KalmanEnergyErr    = 0.;
    m_TkrCal_radlen      = 0.;

    TObjArray::Clear();
}

void TkrTrack::Print(Option_t *option) const {
    TObject::Print(option);
}
/*
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
*/

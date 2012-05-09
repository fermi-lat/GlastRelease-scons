#include "reconRootData/TkrTrack.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

#include "ReconObjectManager.h"

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
    m_rangeEnergy        = 0.;

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

//======================================================
// For Unit Tests
//======================================================

void TkrTrack::Fake( Int_t ievent, UInt_t rank, Float_t randNum ) {

        Float_t f = Float_t(ievent) ;

        UInt_t xgaps = rank * rank ;
        UInt_t ygaps = rank+rank;
        UInt_t x1st = ievent+rank;
        UInt_t y1st = ievent+ievent;
        setNumXGaps(xgaps);
        setNumYGaps(ygaps);
        setNumXFirstGaps(x1st);
        setNumYFirstGaps(y1st);
        Double_t chiSq = f;
        Double_t chiSqSmooth = f*randNum;
        Double_t rms = randNum*randNum;
        Double_t qual = f*f;
        Double_t e = f*2.0;
        Double_t ms = f*6.0;
        setChiSquareFilter(chiSq);
        setChiSquareSmooth(chiSqSmooth);
        setScatter(rms);
        setQuality(qual);
        setKalEnergy(e);
        setKalThetaMS(ms);
}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrTrack::CompareInRange( const TkrTrack & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(NumXGaps) && result ;
    result = COMPARE_IN_RANGE(NumYGaps) && result ;
    result = COMPARE_IN_RANGE(NumXFirstGaps) && result ;
    result = COMPARE_IN_RANGE(NumYFirstGaps) && result ;
    result = COMPARE_IN_RANGE(ChiSquareFilter) && result ;
    result = COMPARE_IN_RANGE(ChiSquareSmooth) && result ;
    result = COMPARE_IN_RANGE(Scatter) && result ;
    result = COMPARE_IN_RANGE(Quality) && result ;
    result = COMPARE_IN_RANGE(KalEnergy) && result ;
    result = COMPARE_IN_RANGE(KalThetaMS) && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

void* TkrTrack::operator new(size_t /*size*/)
{
    TkrTrack* temp = ReconObjectManager::getPointer()->getNewTkrTrack();

    // Since we recycle, make sure these member functions are cleared
    //temp->m_hitCol.clear();

    return temp;
}

void* TkrTrack::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrTrack::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}


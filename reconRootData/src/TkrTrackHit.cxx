#include "reconRootData/TkrTrackHit.h"
#include <iostream>

ClassImp(TkrTrackHit)

TkrTrackHit::TkrTrackHit() 
{
    Clear();
}

void TkrTrackHit::initializeInfo(TkrCluster*           cluster,
                                 commonRootData::TkrId tkrID,
                                 const Double_t        hitZ,
                                 const Double_t        hitEnergy,
                                 const Double_t        hitRadLen,
                                 const Double_t        hitActDist,
                                 const Double_t        hitChiFilter,
                                 const Double_t        hitChiSmooth)
{
    // First make sure we are "clear"
    Clear();

    // Now initialize
    m_cluster          = cluster; 
    m_hitID            = tkrID;
    m_zPlane           = hitZ;
    m_energy           = hitEnergy; 
    m_radLen           = hitRadLen;
    m_activeDist       = hitActDist; 
    m_chiSquareFilter  = hitChiFilter; 
    m_chiSquareRevFit  = 0.; 
    m_chiSquareSmooth  = hitChiSmooth; 

    if (m_cluster != 0)  m_statusBits  = HITISSSD | HITONFIT;
    if (tkrID.hasTray()) m_statusBits |= HASVALIDTKR;
}

void TkrTrackHit::initializeHits(const TkrTrackParams& meas, 
                                 const TkrTrackParams& pred, 
                                 const TkrTrackParams& fit,
                                 const TkrTrackParams& smooth, 
                                 const TkrTrackParams& material)
{
    m_hitMeas   = meas;
    m_hitPred   = pred;
    m_hitFit    = fit;
    m_hitSmooth = smooth;
    m_Qmaterial = material;
}

void TkrTrackHit::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_statusBits       = 0;                       // See StatusBits enumeration above for definitions
    m_cluster          = 0;                       // Pointer to the cluster associated with this hit
    m_hitID            = commonRootData::TkrId(); // Complete TkrId identifying the details of this hit/plane
    m_zPlane           = 0.;                      // Z location of plane
    m_energy           = 0.;                      // Energy of track at this plane
    m_radLen           = 0.;                      // Radiation Lengths encountered from the previous hit
    m_activeDist       = 0.;                      // The distance inside (positive) hit SSD (neg. if outside)
    m_chiSquareFilter  = 0.;                      // hit chi-square at filter stage of fit
    m_chiSquareRevFit  = 0.;                      // hit chi-square at filter stage of fit
    m_chiSquareSmooth  = 0.;                      // hit chi-square at smooth stage of fit

    m_hitMeas          = TkrTrackParams();
    m_hitPred          = TkrTrackParams();
    m_hitFit           = TkrTrackParams();
    m_hitSmooth        = TkrTrackParams();
    m_Qmaterial        = TkrTrackParams();
}

void TkrTrackHit::Print(Option_t *option) const {
    TObject::Print(option);
}


const TVector3 TkrTrackHit::getPoint(TkrTrackHit::ParamType type) const
{
    const TkrTrackParams& hit = getTrackParams(type);

    return TVector3(hit.getxPosition(), hit.getyPosition(),m_zPlane);
}

TVector3 TkrTrackHit::getPoint(TkrTrackHit::ParamType type) 
{
    const TkrTrackParams& hit = getTrackParams(type);

    return TVector3(hit.getxPosition(), hit.getyPosition(),m_zPlane);
}

const TVector3 TkrTrackHit::getDirection(TkrTrackHit::ParamType type) const
{
    const TkrTrackParams& hit = getTrackParams(type);

    TVector3 dir;

    // This assumes that track directions are slopes
	if(m_statusBits & UPWARDS)
    {
		dir = TVector3(hit.getxSlope(), hit.getySlope(), 1.);
	}
	else 
    { 
		dir = TVector3(-hit.getxSlope(), -hit.getySlope(), -1.);
	}

    dir.SetMag(1.);

    // This assumes that track directions are slopes
    return dir;
}

TVector3 TkrTrackHit::getDirection(TkrTrackHit::ParamType type)
{
    const TkrTrackParams& hit = getTrackParams(type);

    TVector3 dir;

    // This assumes that track directions are slopes
	if(m_statusBits & UPWARDS)
    {
		dir = TVector3(hit.getxSlope(), hit.getySlope(), 1.);
	}
	else 
    { 
		dir = TVector3(-hit.getxSlope(), -hit.getySlope(), -1.);
	}

    dir.SetMag(1.);

    return dir;
}

const TkrTrackParams& TkrTrackHit::getTrackParams(TkrTrackHit::ParamType type) const
{
    //
    // Access via this method presumes the track parameters are set. So, REQUIRE the
    // status bits to be valid before returning the hit information
    //
    switch(type)
    {
        case MEASURED:  
        {
            if (!validMeasuredHit()) throw std::invalid_argument("Invalid Measured TkrTrackParams requested");
            return m_hitMeas;
        }
        case PREDICTED: 
        {
            if (!validPredictedHit()) throw std::invalid_argument("Invalid Predicted TkrTrackParams requested");
            return m_hitPred;
        }
        case FILTERED:  
        {
            if (!validFilteredHit()) throw std::invalid_argument("Invalid Filtered TkrTrackParams requested");
            return m_hitFit;
        }
        case REVFIT:  
        {
            if (!validFilteredHit()) throw std::invalid_argument("Invalid Filtered TkrTrackParams requested");
            return m_hitRevFit;
        }
        case SMOOTHED:  
        {
            if (!validSmoothedHit()) throw std::invalid_argument("Invalid Smoothed TkrTrackParams requested");
            return m_hitSmooth;
        }
        case QMATERIAL: 
        {
            if (!validMaterial()) throw std::invalid_argument("Invalid Material TkrTrackParams requested");
            return m_Qmaterial;
        }
        default:
        {
            throw std::invalid_argument("Invalid type of TkrTrackParams hit requested");
        }
    }

    return m_hitMeas;
}

TkrTrackParams& TkrTrackHit::getTrackParams(TkrTrackHit::ParamType type)
{
    //
    // Access to track parameters via this method does not require valid status bit 
    // (since this might be need to set them)
    // 
    switch(type)
    {
        case MEASURED:  return m_hitMeas;
        case PREDICTED: return m_hitPred;
        case FILTERED:  return m_hitFit;
        case REVFIT:    return m_hitRevFit;
        case SMOOTHED:  return m_hitSmooth;
        case QMATERIAL: return m_Qmaterial;
        case UNKNOWN:   
            throw std::invalid_argument("Type UNKNOWN of TkrTrackParams hit requested");

        // If here then something is wrong
        throw std::invalid_argument("Invalid type of TkrTrackParams hit requested");
    }

    return m_hitMeas;
}

const int TkrTrackHit::getParamIndex(TkrTrackHit::SSDDirection meas, TkrTrackParams::ParamType type) const
{
    int posIdx = TkrTrackParams::xPosIdx;
    int slpIdx = TkrTrackParams::xSlpIdx;

    if (   ( meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == commonRootData::TkrId::eMeasureY) 
        || (!meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == commonRootData::TkrId::eMeasureX) )
    {
        posIdx = TkrTrackParams::yPosIdx;
        slpIdx = TkrTrackParams::ySlpIdx;
    }

    return type == TkrTrackParams::Slope ? slpIdx : posIdx;
}

int TkrTrackHit::getParamIndex(TkrTrackHit::SSDDirection meas, TkrTrackParams::ParamType type)
{
    int posIdx = TkrTrackParams::xPosIdx;
    int slpIdx = TkrTrackParams::xSlpIdx;

    if (   ( meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == commonRootData::TkrId::eMeasureY) 
        || (!meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == commonRootData::TkrId::eMeasureX) )
    {
        posIdx = TkrTrackParams::yPosIdx;
        slpIdx = TkrTrackParams::ySlpIdx;
    }

    return type == TkrTrackParams::Slope ? slpIdx : posIdx;
}

const double TkrTrackHit::getMeasuredPosition(TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDMEASURED, TkrTrackParams::Position));
}

const double TkrTrackHit::getMeasuredSlope(TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDMEASURED, TkrTrackParams::Slope));
}

const double TkrTrackHit::getNonMeasuredPosition(TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDNONMEASURED, TkrTrackParams::Position));
}

const double TkrTrackHit::getNonMeasuredSlope(TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDNONMEASURED, TkrTrackParams::Slope));
}

double TkrTrackHit::getCoordinate(const TkrTrackParams& params, int coord) 
              {return params(coord);}

const double TkrTrackHit::getCoordinate(const TkrTrackParams& params, int coord) const 
              {return params(coord);}


//======================================================
// For Unit Tests
//======================================================

void TkrTrackHit::Fake( Int_t ievent, UInt_t rank, Float_t randNum ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrTrackHit::CompareInRange( const TkrTrackHit & ref, const std::string & name ) const {

    bool result = true ;

    //result = COMPARE_IN_RANGE(Energy) && result ;

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


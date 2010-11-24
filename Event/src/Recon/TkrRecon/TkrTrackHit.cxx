
#include "Event/Recon/TkrRecon/TkrTrackHit.h"
#include <stdexcept>


void Event::TkrTrackHit::clean()
{
    Event::TkrTrackParams measParams = m_tkrParams[MEASURED];

    m_tkrParams.clear();

    m_tkrParams[MEASURED] = measParams;

    m_statusBits &= HITONFIT | HASMEASURED | HASVALIDTKR;

    m_chiSquareFilter = 0.;
    m_chiSquareSmooth = 0.;
}

void Event::TkrTrackHit::clear()
{
    m_tkrParams.clear();

    m_statusBits &= HITONFIT | HASVALIDTKR;

    m_chiSquareFilter = 0.;
    m_chiSquareSmooth = 0.;
}

const Point Event::TkrTrackHit::getPoint(TkrTrackHit::ParamType type) const
{
    const TkrTrackParams& hit = getTrackParams(type);

    return Point(hit.getxPosition(), hit.getyPosition(),m_zPlane);
}

Point Event::TkrTrackHit::getPoint(TkrTrackHit::ParamType type) 
{
    const TkrTrackParams& hit = getTrackParams(type);

    return Point(hit.getxPosition(), hit.getyPosition(),m_zPlane);
}

const Vector Event::TkrTrackHit::getDirection(TkrTrackHit::ParamType type) const
{
    TkrTrackParams hit = getTrackParams(type);

    if (m_statusBits & HITHASKINKANG)
    {
        int    measIdx = getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Position);
        double measAng = atan(hit(measIdx));

        measAng += m_kinkAngle;

        hit(measIdx) = tan(measAng);
    }

    // This assumes that track directions are slopes
    return Vector(-hit.getxSlope(), -hit.getySlope(), -1.).unit();
}

Vector Event::TkrTrackHit::getDirection(TkrTrackHit::ParamType type)
{
    TkrTrackParams hit = getTrackParams(type);

    if (m_statusBits & HITHASKINKANG)
    {
        int    measIdx = getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Position);
        double measAng = atan(hit(measIdx));

        measAng += m_kinkAngle;

        hit(measIdx) = tan(measAng);
    }

    // This assumes that track directions are slopes
	if(m_statusBits & UPWARDS){
		return Vector(hit.getxSlope(), hit.getySlope(), 1.).unit();
	}
	else { 
		return 	Vector(-hit.getxSlope(), -hit.getySlope(), -1.).unit();
	}
}

const Event::TkrTrackParams& Event::TkrTrackHit::getTrackParams(TkrTrackHit::ParamType type) const
{
    //
    // Access via this method presumes the track parameters are set. So, REQUIRE the
    // status bits to be valid before returning the hit information
    //
    TkrParamsMap::const_iterator mapIter = m_tkrParams.find(type);

    if (mapIter == m_tkrParams.end()) 
    {
        throw std::invalid_argument("Invalid Measured TkrTrackParams requested");
    }

    return mapIter->second;
}

Event::TkrTrackParams& Event::TkrTrackHit::getTrackParams(TkrTrackHit::ParamType type)
{
    //
    // Access to track parameters via this method does not require valid status bit 
    // (since this might be need to set them)
    // 
    TkrParamsMap::iterator mapIter = m_tkrParams.find(type);

    if (mapIter == m_tkrParams.end()) 
    {
        // They don't yet exist so create them
        Event::TkrTrackParams params = Event::TkrTrackParams();

        m_tkrParams[type] = params;

        mapIter = m_tkrParams.find(type);
    }

    return mapIter->second;
}

void Event::TkrTrackHit::setTrackParams(ITkrTrackParamsAccess& access, TkrTrackHit::ParamType type) 
{
    // Don't forget to update the status bits...
    unsigned int statusBits = 0;

    Event::TkrTrackParams& params = getTrackParams(type);

    access.setParams(&params);

    m_tkrParams[type] = params;

    statusBits = 1 << (type + 1);

/*
    // Switch on the type of hit to set
    switch(type)
    {
        case MEASURED:
        {
            access.setParams(&m_hitMeas);
            statusBits = HASMEASURED;
            break;
        }
        case PREDICTED: 
        {
            access.setParams(&m_hitPred);
            statusBits = HASPREDICTED;
            break;
        }
        case FILTERED:  
        {
            access.setParams(&m_hitFit);
            statusBits = HASFILTERED;
            break;
        }
        case REVFIT:  
        {
            access.setParams(&m_hitRevFit);
            statusBits = HASFILTERED;
            break;
        }
        case SMOOTHED:
        {
            access.setParams(&m_hitSmooth);
            statusBits = HASSMOOTHED;
            break;
        }
        case QMATERIAL: 
        {
            access.setParams(&m_Qmaterial);
            statusBits = HASMATERIAL;
        }
        case UNKNOWN:
        default: {}
    }
*/

    m_statusBits |= statusBits;

    return;
}

const int Event::TkrTrackHit::getParamIndex(TkrTrackHit::SSDDirection meas, TkrTrackParams::ParamType type) const
{
    int posIdx = Event::TkrTrackParams::xPosIdx;
    int slpIdx = Event::TkrTrackParams::xSlpIdx;

    if (   ( meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == idents::TkrId::eMeasureY) 
        || (!meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == idents::TkrId::eMeasureX) )
    {
        posIdx = Event::TkrTrackParams::yPosIdx;
        slpIdx = Event::TkrTrackParams::ySlpIdx;
    }

    return type == TkrTrackParams::Slope ? slpIdx : posIdx;
}

int Event::TkrTrackHit::getParamIndex(TkrTrackHit::SSDDirection meas, TkrTrackParams::ParamType type)
{
    int posIdx = Event::TkrTrackParams::xPosIdx;
    int slpIdx = Event::TkrTrackParams::xSlpIdx;

    if (   ( meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == idents::TkrId::eMeasureY) 
        || (!meas == TkrTrackHit::SSDMEASURED && m_hitID.getView() == idents::TkrId::eMeasureX) )
    {
        posIdx = Event::TkrTrackParams::yPosIdx;
        slpIdx = Event::TkrTrackParams::ySlpIdx;
    }

    return type == TkrTrackParams::Slope ? slpIdx : posIdx;
}

const double Event::TkrTrackHit::getMeasuredPosition(Event::TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDMEASURED, TkrTrackParams::Position));
}

const double Event::TkrTrackHit::getMeasuredSlope(Event::TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDMEASURED, TkrTrackParams::Slope));
}

const double Event::TkrTrackHit::getNonMeasuredPosition(Event::TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDNONMEASURED, TkrTrackParams::Position));
}

const double Event::TkrTrackHit::getNonMeasuredSlope(Event::TkrTrackHit::ParamType type) const
{
    return getCoordinate(getTrackParams(type), getParamIndex(SSDNONMEASURED, TkrTrackParams::Slope));
}

std::ostream& Event::TkrTrackHit::fillStream( std::ostream& s ) const 
{ 
  s << "StatusBit   : "<< getStatusBits()<<"\n"
    << "Z Plane     : "<< getZPlane()<<"\n"
    <<" Energy      : "<< getEnergy() <<"\n"
    <<" Energy      : "<< getRadLen()
    <<" Active Dist : "<< getActiveDist()
    <<" Filter Chi2 : "<< getChiSquareFilter()
    <<" Smooth Chi2 : "<< getChiSquareSmooth();
  return s; 
}

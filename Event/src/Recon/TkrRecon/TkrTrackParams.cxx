
#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include <stdexcept>

using namespace Event;

Event::TkrTrackParams::TkrTrackParams()
{
    initDataMembers();

    return;
}

Event::TkrTrackParams::TkrTrackParams(const TkrTrackParams& right)
{
    *this = right;

    return;
}

Event::TkrTrackParams::TkrTrackParams(ITkrTrackParamsAccess& access)
{
    access.setParams(this);

    return;
}


//Event::TkrTrackParams::TkrTrackParams(ParamType t,
//                   double xPosition, double xSlope, double yPosition, double ySlope,
//                   double xPosxPos, double xPosxSlp, double xPosyPos, double xPosySlp,
//                   double xSlpxSlp, double xSlpyPos, double xSlpySlp,
//                   double yPosyPos, double yPosySlp,
//                   double ySlpySlp);

//Event::TkrTrackParams::TkrTrackParams (const TkrTrackParams& right);

void Event::TkrTrackParams::initDataMembers()
{
    /// Track parameters
    m_xPosition  = 0.;
    m_xSlope = 0.;
    m_yPosition  = 0.;
    m_ySlope = 0.;

    /// Track parameter error matrix elements
    m_xPos_xPos  = 0.;     // Cov(1,1) = xPositionErr * xPositionErr
    m_xPos_xSlp  = 0.;     // Cov(1,2) = Cov(2,1) = xPositionErr * xSlopeErr
    m_xPos_yPos  = 0.;     // Cov(1,3) = Cov(3,1) = xPositionErr * yPositionErr
    m_xPos_ySlp  = 0.;     // Cov(1,4) = Cov(4,1) = xPositionERr * ySlopeErr
    m_xSlp_xSlp  = 0.;     // Cov(2,2) = xSlopeErr * xSlopeErr
    m_xSlp_yPos  = 0.;     // Cov(2,3) = Cov(3,2) = xSlopeErr * yPositionErr
    m_xSlp_ySlp  = 0.;     // Cov(2,4) = Cov(4,2) = xSlopeErr * ySlopeErr
    m_yPos_yPos  = 0.;     // Cov(3,3) = yPositionErr * yPositionErr
    m_yPos_ySlp  = 0.;     // Cov(3,4) = Cov(4,3) = yPositionErr * ySlopeErr
    m_ySlp_ySlp  = 0.;     // Cov(4,4) = ySlopeErr * ySlopeErr

    return;
}

double& Event::TkrTrackParams::operator()(const int &i)
{
    if (i < xPosIdx || i > ySlpIdx) 
          throw std::invalid_argument("Invalid index for TkrTrackParams");

    switch(i)
    {
        case 1: return m_xPosition;
        case 2: return m_xSlope;
        case 3: return m_yPosition;
    }

    return m_ySlope;
}

const double& Event::TkrTrackParams::operator()(const int &i) const
{
    if (i < xPosIdx || i > ySlpIdx) 
          throw std::invalid_argument("Invalid index for TkrTrackParams");

    switch(i)
    {
        case 1: return m_xPosition;
        case 2: return m_xSlope;
        case 3: return m_yPosition;
    }

    return m_ySlope;
}

double& Event::TkrTrackParams::operator()(const int &i, const int &j)
{
    if (i < xPosIdx || j < xPosIdx || i > ySlpIdx || j > ySlpIdx) 
          throw std::invalid_argument("Invalid index for TkrTrackParams");

    // yuk...
    if       (i == 1 && j == 1)                        return m_xPos_xPos;
    else if ((i == 1 && j == 2) || (i == 2 && j == 1)) return m_xPos_xSlp;
    else if ((i == 1 && j == 3) || (i == 3 && j == 1)) return m_xPos_yPos;
    else if ((i == 1 && j == 4) || (i == 4 && j == 1)) return m_xPos_ySlp;
    else if  (i == 2 && j == 2)                        return m_xSlp_xSlp;
    else if ((i == 2 && j == 3) || (i == 3 && j == 2)) return m_xSlp_yPos;
    else if ((i == 2 && j == 4) || (i == 4 && j == 2)) return m_xSlp_ySlp;
    else if  (i == 3 && j == 3)                        return m_yPos_yPos;
    else if ((i == 3 && j == 4) || (i == 4 && j == 3)) return m_yPos_ySlp;

    return m_ySlp_ySlp;
}

const double& Event::TkrTrackParams::operator()(const int &i, const int &j) const
{
    if (i < xPosIdx || j < xPosIdx || i > ySlpIdx || j > ySlpIdx) 
          throw std::invalid_argument("Invalid index for TkrTrackParams");

    // yuk...
    if       (i == 1 && j == 1)                        return m_xPos_xPos;
    else if ((i == 1 && j == 2) || (i == 2 && j == 1)) return m_xPos_xSlp;
    else if ((i == 1 && j == 3) || (i == 3 && j == 1)) return m_xPos_yPos;
    else if ((i == 1 && j == 4) || (i == 4 && j == 1)) return m_xPos_ySlp;
    else if  (i == 2 && j == 2)                        return m_xSlp_xSlp;
    else if ((i == 2 && j == 3) || (i == 3 && j == 2)) return m_xSlp_yPos;
    else if ((i == 2 && j == 4) || (i == 4 && j == 2)) return m_xSlp_ySlp;
    else if  (i == 3 && j == 3)                        return m_yPos_yPos;
    else if ((i == 3 && j == 4) || (i == 4 && j == 3)) return m_yPos_ySlp;

    return m_ySlp_ySlp;
}


std::ostream& Event::TkrTrackParams::fillStream( std::ostream& s ) const 
{ 
  s << getxPosition() <<" "<<getxSlope()<<" "<<getyPosition() <<" "<<getySlope()<<"\n"
    << getxPosxPos() <<" "<<getxPosxSlp()<<" "<<getxPosyPos()<<" "<<getxPosySlp()<<"\n"
    <<"          "<<getxSlpxSlp()<<" "<<getxSlpyPos()<<" "<<getxSlpySlp()<<"\n"
    <<"          "<<"          "<<getyPosyPos()<<" "<<getyPosySlp()<<"\n"
    <<"          "<<"          "<<"          "<<getySlpySlp();
  return s; 
}



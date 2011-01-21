/** @file TkrTrackParams.h
*
* @author Bill Atwood, Leon Rochester, Johann Cohen, Tracy Usher
*
* $Header$

*/
#ifndef TkrTrackParams_H
#define TkrTrackParams_H


#include <iostream>

namespace Event { //Namespace Event

/** 
* @class TkrTrackParams
*
* @brief Gaudi TDS class to store the track parameters and associated covariance
*        matrix. This class is purely data storage and its data is accessed through
*        methods provided in the TkrFitPlane class which owns it.
* 
*        This updates the old TkrFitHit class.
*/

// An interface class to be used for filling and retrieving data by "special" helper classes
class TkrTrackParams;

class ITkrTrackParamsAccess
{
public:
    virtual void setParams(TkrTrackParams* params) = 0;
    virtual void getParams(TkrTrackParams* params) = 0;
};

class TkrTrackParams
{
public:
    /// Default constructor
    TkrTrackParams();

    /// Construct from an external object
    TkrTrackParams(ITkrTrackParamsAccess& access);

    /// Direct construction from all the elements (the old fashioned way)
    TkrTrackParams(double xPosition, double xSlope, double yPosition, double ySlope,
                   double xPosxPos, double xPosxSlp, double xPosyPos, double xPosySlp,
                   double xSlpxSlp, double xSlpyPos, double xSlpySlp,
                   double yPosyPos, double yPosySlp,
                   double ySlpySlp);

    /// Copy constructor
    TkrTrackParams (const TkrTrackParams& right);

    /// Direct access to parameters
    inline const double getxPosition()  const {return m_xPosition; }
    inline const double getxSlope()     const {return m_xSlope;    }
    inline const double getyPosition()  const {return m_yPosition; }
    inline const double getySlope()     const {return m_ySlope;    }

    /// Direct access to errors
    inline const double getxPosxPos()   const {return m_xPos_xPos; }
    inline const double getxPosxSlp()   const {return m_xPos_xSlp; }
    inline const double getxPosyPos()   const {return m_xPos_yPos; }
    inline const double getxPosySlp()   const {return m_xPos_ySlp; }
    inline const double getxSlpxSlp()   const {return m_xSlp_xSlp; }
    inline const double getxSlpyPos()   const {return m_xSlp_yPos; }
    inline const double getxSlpySlp()   const {return m_xSlp_ySlp; }
    inline const double getyPosyPos()   const {return m_yPos_yPos; }
    inline const double getyPosySlp()   const {return m_yPos_ySlp; }
    inline const double getySlpySlp()   const {return m_ySlp_ySlp; }

    /// Define track parameter types and indices for look up with the array operators
    enum ParamType  {Position, Slope};
    enum ParamIndex {xPosIdx = 1, xSlpIdx = 2, yPosIdx = 3, ySlpIdx = 4};

    /// Define an ( ) operator (allows read/write - indexing from 1!!)
    double& operator()(const int &i);
    const double& operator()(const int &i) const;
    double& operator()(const int &i, const int &j);
    const double& operator()(const int &i, const int &j) const;

    /// Set parameters
    inline void setxPosition(const double& val)  {m_xPosition = val; }
    inline void setxSlope(const double& val)     {m_xSlope = val;    }
    inline void setyPosition(const double& val)  {m_yPosition = val; }
    inline void setySlope(const double& val)     {m_ySlope = val;    }

    /// Set errors
    inline void setxPosxPos(const double& val)   {m_xPos_xPos = val; }
    inline void setxPosxSlp(const double& val)   {m_xPos_xSlp = val; }
    inline void setxPosyPos(const double& val)   {m_xPos_yPos = val; }
    inline void setxPosySlp(const double& val)   {m_xPos_ySlp = val; }
    inline void setxSlpxSlp(const double& val)   {m_xSlp_xSlp = val; }
    inline void setxSlpyPos(const double& val)   {m_xSlp_yPos = val; }
    inline void setxSlpySlp(const double& val)   {m_xSlp_ySlp = val; }
    inline void setyPosyPos(const double& val)   {m_yPos_yPos = val; }
    inline void setyPosySlp(const double& val)   {m_yPos_ySlp = val; }
    inline void setySlpySlp(const double& val)   {m_ySlp_ySlp = val; }

    std::ostream& fillStream( std::ostream& s ) const;
    friend std::ostream& operator<< ( std::ostream& s, const TkrTrackParams& obj ) 
      {
        return obj.fillStream(s);
      }

private:
    /// Private intialization method
    void initDataMembers();
  
    /// Track parameters
    double m_xPosition;     // x position parameter
    double m_xSlope;        // x "Slope" (x slope or tan(theta x)
    double m_yPosition;     // y position
    double m_ySlope;        // y "Slope" (y slope or tan(theta y)

    /// Track parameter error matrix elements
    double m_xPos_xPos;     // Cov(1,1) = xPositionErr * xPositionErr
    double m_xPos_xSlp;     // Cov(1,2) = Cov(2,1) = xPositionErr * xSlopeErr
    double m_xPos_yPos;     // Cov(1,3) = Cov(3,1) = xPositionErr * yPositionErr
    double m_xPos_ySlp;     // Cov(1,4) = Cov(4,1) = xPositionERr * ySlopeErr
    double m_xSlp_xSlp;     // Cov(2,2) = xSlopeErr * xSlopeErr
    double m_xSlp_yPos;     // Cov(2,3) = Cov(3,2) = xSlopeErr * yPositionErr
    double m_xSlp_ySlp;     // Cov(2,4) = Cov(4,2) = xSlopeErr * ySlopeErr
    double m_yPos_yPos;     // Cov(3,3) = yPositionErr * yPositionErr
    double m_yPos_ySlp;     // Cov(3,4) = Cov(4,3) = yPositionErr * ySlopeErr
    double m_ySlp_ySlp;     // Cov(4,4) = ySlopeErr * ySlopeErr
};

}; //Namespace Event


#endif


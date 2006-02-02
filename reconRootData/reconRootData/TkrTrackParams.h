#ifndef ROOT_TKRTRACKPARAMS_H
#define ROOT_TKRTRACKPARAMS_H

#include "TObject.h"
/** 
 * @class TkrTrackParams
 * @brief Root object containing fit track parameters for a given plane
 * 
 * This class should be a duplicate of Event::TkrTrackParams
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

class TkrTrackParams : public TObject
{
public:

    TkrTrackParams();

    TkrTrackParams(Double32_t xPosition, Double32_t xSlope, Double32_t yPosition, Double32_t ySlope,
                   Double32_t xPosxPos, Double32_t xPosxSlp, Double32_t xPosyPos, Double32_t xPosySlp,
                   Double32_t xSlpxSlp, Double32_t xSlpyPos, Double32_t xSlpySlp,
                   Double32_t yPosyPos, Double32_t yPosySlp,
                   Double32_t ySlpySlp);

    TkrTrackParams(const TkrTrackParams& params);

    virtual ~TkrTrackParams() {};

    /// resets all member variables to default values
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrTrackParams &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    TkrTrackParams& operator=(const TkrTrackParams& old);

    /// Provide an initialization method (in case of default constructor)
    void initialize(Double32_t xPosition, Double32_t xSlope, Double32_t yPosition, Double32_t ySlope,
                    Double32_t xPosxPos, Double32_t xPosxSlp, Double32_t xPosyPos, Double32_t xPosySlp,
                    Double32_t xSlpxSlp, Double32_t xSlpyPos, Double32_t xSlpySlp,
                    Double32_t yPosyPos, Double32_t yPosySlp,
                    Double32_t ySlpySlp);
    void initialize(const TkrTrackParams& params);

    /// Direct access to parameters
    inline const Double32_t getxPosition()  const {return m_xPosition; }
    inline const Double32_t getxSlope()     const {return m_xSlope;    }
    inline const Double32_t getyPosition()  const {return m_yPosition; }
    inline const Double32_t getySlope()     const {return m_ySlope;    }

    /// Direct access to errors
    inline const Double32_t getxPosxPos()   const {return m_xPos_xPos; }
    inline const Double32_t getxPosxSlp()   const {return m_xPos_xSlp; }
    inline const Double32_t getxPosyPos()   const {return m_xPos_yPos; }
    inline const Double32_t getxPosySlp()   const {return m_xPos_ySlp; }
    inline const Double32_t getxSlpxSlp()   const {return m_xSlp_xSlp; }
    inline const Double32_t getxSlpyPos()   const {return m_xSlp_yPos; }
    inline const Double32_t getxSlpySlp()   const {return m_xSlp_ySlp; }
    inline const Double32_t getyPosyPos()   const {return m_yPos_yPos; }
    inline const Double32_t getyPosySlp()   const {return m_yPos_ySlp; }
    inline const Double32_t getySlpySlp()   const {return m_ySlp_ySlp; }

    /// Define track parameter types and indices for look up with the array operators
    enum ParamType  {Position, Slope};
    enum ParamIndex {xPosIdx = 1, xSlpIdx = 2, yPosIdx = 3, ySlpIdx = 4};

    /// Define an ( ) operator (allows read/write - indexing from 1!!)
    Double32_t&       operator()(const int &i);
    const Double32_t& operator()(const int &i) const;
    Double32_t&       operator()(const int &i, const int &j);
    const Double32_t& operator()(const int &i, const int &j) const;

private:
  
    /// Track parameters
    Double32_t m_xPosition;     // x position parameter
    Double32_t m_xSlope;        // x "Slope" (x slope or tan(theta x)
    Double32_t m_yPosition;     // y position
    Double32_t m_ySlope;        // y "Slope" (y slope or tan(theta y)

    /// Track parameter error matrix elements
    Double32_t m_xPos_xPos;     // Cov(1,1) = xPositionErr * xPositionErr
    Double32_t m_xPos_xSlp;     // Cov(1,2) = Cov(2,1) = xPositionErr * xSlopeErr
    Double32_t m_xPos_yPos;     // Cov(1,3) = Cov(3,1) = xPositionErr * yPositionErr
    Double32_t m_xPos_ySlp;     // Cov(1,4) = Cov(4,1) = xPositionERr * ySlopeErr
    Double32_t m_xSlp_xSlp;     // Cov(2,2) = xSlopeErr * xSlopeErr
    Double32_t m_xSlp_yPos;     // Cov(2,3) = Cov(3,2) = xSlopeErr * yPositionErr
    Double32_t m_xSlp_ySlp;     // Cov(2,4) = Cov(4,2) = xSlopeErr * ySlopeErr
    Double32_t m_yPos_yPos;     // Cov(3,3) = yPositionErr * yPositionErr
    Double32_t m_yPos_ySlp;     // Cov(3,4) = Cov(4,3) = yPositionErr * ySlopeErr
    Double32_t m_ySlp_ySlp;     // Cov(4,4) = ySlopeErr * ySlopeErr

    ClassDef(TkrTrackParams,2)

};
#endif

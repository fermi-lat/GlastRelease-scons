#include "reconRootData/TkrTrackParams.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrTrackParams)

TkrTrackParams::TkrTrackParams() 
{
    Clear();
}


TkrTrackParams::TkrTrackParams(Double_t xPosition, Double_t xSlope, Double_t yPosition, Double_t ySlope,
                   Double_t xPosxPos, Double_t xPosxSlp, Double_t xPosyPos, Double_t xPosySlp,
                   Double_t xSlpxSlp, Double_t xSlpyPos, Double_t xSlpySlp,
                   Double_t yPosyPos, Double_t yPosySlp,
                   Double_t ySlpySlp) 
{
    initialize(xPosition, xSlope, yPosition, ySlope,
               xPosxPos, xPosxSlp, xPosyPos, xPosySlp,
               xSlpxSlp, xSlpyPos, xSlpySlp,
               yPosyPos, yPosySlp,
               ySlpySlp);
}

TkrTrackParams::TkrTrackParams(const TkrTrackParams& params) : TObject(params)
{
    initialize(params);
}

void TkrTrackParams::Clear(Option_t* /* option */) 
{
   /// Track parameters
    m_xPosition  = 0.;
    m_xSlope     = 0.;
    m_yPosition  = 0.;
    m_ySlope     = 0.;

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
}

void TkrTrackParams::Print(Option_t *option) const 
{
    TObject::Print(option);
    using namespace std;

    cout << "(xpos, xslope) = (" << m_xPosition << ", " << m_xSlope << ")" << endl;
    cout << "(ypos, yslope) = (" << m_yPosition << ", " << m_ySlope << ")" << endl;

}

TkrTrackParams& TkrTrackParams::operator=(const TkrTrackParams& params)
{
    initialize(params);

    return *this;
}

void TkrTrackParams::initialize(Double_t xPosition, Double_t xSlope, Double_t yPosition, Double_t ySlope,
                    Double_t xPosxPos, Double_t xPosxSlp, Double_t xPosyPos, Double_t xPosySlp,
                    Double_t xSlpxSlp, Double_t xSlpyPos, Double_t xSlpySlp,
                    Double_t yPosyPos, Double_t yPosySlp,
                    Double_t ySlpySlp)
{
   /// Track parameters
    m_xPosition  = xPosition;
    m_xSlope     = xSlope;
    m_yPosition  = yPosition;
    m_ySlope     = ySlope;

    /// Track parameter error matrix elements
    m_xPos_xPos  = xPosxPos;     // Cov(1,1) = xPositionErr * xPositionErr
    m_xPos_xSlp  = xPosxSlp;     // Cov(1,2) = Cov(2,1) = xPositionErr * xSlopeErr
    m_xPos_yPos  = xPosyPos;     // Cov(1,3) = Cov(3,1) = xPositionErr * yPositionErr
    m_xPos_ySlp  = xPosySlp;     // Cov(1,4) = Cov(4,1) = xPositionERr * ySlopeErr
    m_xSlp_xSlp  = xSlpxSlp;     // Cov(2,2) = xSlopeErr * xSlopeErr
    m_xSlp_yPos  = xSlpyPos;     // Cov(2,3) = Cov(3,2) = xSlopeErr * yPositionErr
    m_xSlp_ySlp  = xSlpySlp;     // Cov(2,4) = Cov(4,2) = xSlopeErr * ySlopeErr
    m_yPos_yPos  = yPosyPos;     // Cov(3,3) = yPositionErr * yPositionErr
    m_yPos_ySlp  = yPosySlp;     // Cov(3,4) = Cov(4,3) = yPositionErr * ySlopeErr
    m_ySlp_ySlp  = ySlpySlp;     // Cov(4,4) = ySlopeErr * ySlopeErr
}  

void TkrTrackParams::initialize(const TkrTrackParams& params)
{
    initialize(params.getxPosition(), params.getxSlope(), params.getyPosition(), params.getySlope(),
               params.getxPosxPos(), params.getxPosxSlp(), params.getxPosyPos(), params.getxPosySlp(),
               params.getxSlpxSlp(), params.getxSlpyPos(), params.getxSlpySlp(), 
               params.getyPosyPos(), params.getyPosySlp(),
               params.getySlpySlp() );
}  


const Double_t& TkrTrackParams::operator() (const int &element) const
{
    switch(element)
    {
        case xPosIdx: return m_xPosition;
        case xSlpIdx: return m_xSlope;
        case yPosIdx: return m_yPosition;
    }

    return m_ySlope;
}

Double_t& TkrTrackParams::operator() (const int &element)
{
    switch(element)
    {
        case xPosIdx: return m_xPosition;
        case xSlpIdx: return m_xSlope;
        case yPosIdx: return m_yPosition;
    }

    return m_ySlope;
}

Double_t& TkrTrackParams::operator()(const int &i, const int &j)
{
    // yuk...
    if       (i == xPosIdx && j == xPosIdx)                                    return m_xPos_xPos;
    else if ((i == xPosIdx && j == xSlpIdx) || (i == xSlpIdx && j == xPosIdx)) return m_xPos_xSlp;
    else if ((i == xPosIdx && j == yPosIdx) || (i == yPosIdx && j == xPosIdx)) return m_xPos_yPos;
    else if ((i == xPosIdx && j == ySlpIdx) || (i == ySlpIdx && j == xPosIdx)) return m_xPos_ySlp;
    else if  (i == xSlpIdx && j == xSlpIdx)                                    return m_xSlp_xSlp;
    else if ((i == xSlpIdx && j == yPosIdx) || (i == yPosIdx && j == xSlpIdx)) return m_xSlp_yPos;
    else if ((i == xSlpIdx && j == ySlpIdx) || (i == ySlpIdx && j == xSlpIdx)) return m_xSlp_ySlp;
    else if  (i == yPosIdx && j == yPosIdx)                                    return m_yPos_yPos;
    else if ((i == yPosIdx && j == ySlpIdx) || (i == ySlpIdx && j == yPosIdx)) return m_yPos_ySlp;

    return m_ySlp_ySlp;
}

const Double_t& TkrTrackParams::operator()(const int &i, const int &j) const
{
    // yuk...
    if       (i == xPosIdx && j == xPosIdx)                                    return m_xPos_xPos;
    else if ((i == xPosIdx && j == xSlpIdx) || (i == xSlpIdx && j == xPosIdx)) return m_xPos_xSlp;
    else if ((i == xPosIdx && j == yPosIdx) || (i == yPosIdx && j == xPosIdx)) return m_xPos_yPos;
    else if ((i == xPosIdx && j == ySlpIdx) || (i == ySlpIdx && j == xPosIdx)) return m_xPos_ySlp;
    else if  (i == xSlpIdx && j == xSlpIdx)                                    return m_xSlp_xSlp;
    else if ((i == xSlpIdx && j == yPosIdx) || (i == yPosIdx && j == xSlpIdx)) return m_xSlp_yPos;
    else if ((i == xSlpIdx && j == ySlpIdx) || (i == ySlpIdx && j == xSlpIdx)) return m_xSlp_ySlp;
    else if  (i == yPosIdx && j == yPosIdx)                                    return m_yPos_yPos;
    else if ((i == yPosIdx && j == ySlpIdx) || (i == ySlpIdx && j == yPosIdx)) return m_yPos_ySlp;

    return m_ySlp_ySlp;
}


//======================================================
// For Unit Tests
//======================================================

void TkrTrackParams::Fake( Int_t ievent, UInt_t /*rank*/, Float_t randNum ) {

    Float_t f = Float_t(ievent) ;

    (*this)(1) = f ;
    (*this)(2) = f*f ;
    (*this)(3) = f*randNum ;
    (*this)(4) = randNum*randNum ;
    
    (*this)(1,1) = f ;
    (*this)(1,2) = f*2.  ;
    (*this)(1,3) = f*3.  ;
    (*this)(1,4) = f*4.  ;
    (*this)(2,2) = f*5.  ;
    (*this)(2,3) = f*6.  ;
    (*this)(2,4) = f*7.  ;
    (*this)(3,3) = f*8.  ;
    (*this)(3,4) = f*9.  ;
    (*this)(4,4) = f*10. ;

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrTrackParams::CompareInRange( const TkrTrackParams & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(xPosition) && result ;
    result = COMPARE_IN_RANGE(xSlope) && result ;
    result = COMPARE_IN_RANGE(yPosition) && result ;
    result = COMPARE_IN_RANGE(ySlope) && result ;

    result = COMPARE_IN_RANGE(xPosxPos) && result ;
    result = COMPARE_IN_RANGE(xPosxSlp) && result ;
    result = COMPARE_IN_RANGE(xPosyPos) && result ;
    result = COMPARE_IN_RANGE(xPosySlp) && result ;
    result = COMPARE_IN_RANGE(xSlpxSlp) && result ;
    result = COMPARE_IN_RANGE(xSlpyPos) && result ;
    result = COMPARE_IN_RANGE(xSlpySlp) && result ;
    result = COMPARE_IN_RANGE(yPosyPos) && result ;
    result = COMPARE_IN_RANGE(yPosySlp) && result ;
    result = COMPARE_IN_RANGE(ySlpySlp) && result ;

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


void* TkrTrackParams::operator new(size_t /*size*/)
{
    TkrTrackParams* temp = ReconObjectManager::getPointer()->getNewTkrTrackParams();

    // Since we recycle, make sure these member functions are cleared
    //temp->m_hitCol.clear();

    return temp;
}

void* TkrTrackParams::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrTrackParams::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}

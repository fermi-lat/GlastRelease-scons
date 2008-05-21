// File and Version Information:
//      $Header$
//
// Description:
//  Implement the ObfFilterStatus class for storing to root the output of the OnboardFilter
//                                                
// Author:  Tracy Usher                    

#include "digiRootData/ObfFilterTrack.h"
#include "Riostream.h"

ClassImp(ObfFilterTrack)


/// Fill the ASCII output stream
void ObfFilterTrack::Print(Option_t *option) const
{
    std::cout << "class ObfFilterTrack"
    << " :"
    << "\n    Number hits X projection  = " << m_nXhits
    << "\n    Number hits Y projection  = " << m_nYhits
    << "\n    Intercept of X projection = " << m_xInt
    << "\n    Intercept of Y projection = " << m_yInt
    << "\n    Z Coordinate at intercept = " << m_z
    << "\n    Slope of X projection     = " << m_slpXZ
    << "\n    Slope of Y projection     = " << m_slpYZ
    << std::endl;

    return;
}

void ObfFilterTrack::initialize(Int_t xHits, Int_t yHits, Float_t xInt, Float_t yInt, Float_t z, Float_t slpXZ, Float_t slpYZ)
{
    m_nXhits = xHits;
    m_nYhits = yHits;
    m_xInt   = xInt;
    m_yInt   = yInt;
    m_z      = z;
    m_slpXZ  = slpXZ;
    m_slpYZ  = slpYZ;
}


void ObfFilterTrack::Clear(Option_t *) 
{
    initialize();
    return;
}

// dummy data, just for tests
void ObfFilterTrack::Fake( Int_t /*ievent*/, UInt_t rank, Float_t randNum ) {

    Clear() ;
/*    
    Float_t f = Float_t(rank) ;
    Float_t fr = f*randNum ;

    ObfFilterStatusPoint* point = new ObfFilterStatusPoint();
    point->Fake(0, rank, randNum);
    addMcPoint(point);

    point = new ObfFilterStatusPoint();
    point->Fake(0, 2*rank, randNum);
    addMcPoint(point);
*/
}

//#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t ObfFilterTrack::CompareInRange( const ObfFilterTrack & ref, const std::string & name ) const {

    Bool_t result = true ;
    
//    result = COMPARE_IN_RANGE(ObfFilterStatusCharge) && result 

//    const McParticle *myPart = getMcParticle();
//    const McParticle *refPart = ref.getMcParticle();
//    result = rootdatautil::CompareInRange(*myPart,*refPart,"McParticle") && result ;


    //result = COMPARE_IN_RANGE(McParticle) && result ;  need to pass reference
/*
    result = rootdatautil::TObjArrayCompareInRange<ObfFilterStatusPoint>(getMcPointCol(),ref.getMcPointCol(),"McPoint") && result ;
      
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
*/
    return result ;

}


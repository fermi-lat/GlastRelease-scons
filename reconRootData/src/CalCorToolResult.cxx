
#include <reconRootData/CalCorToolResult.h>
#include <commonRootData/RootDataUtil.h>
#include <iostream>

ClassImp(CalCorToolResult)

void CalCorToolResult::init
  ( const TString & correctionName,
    UInt_t statusBits,
    Double_t chiSquare,
    const CalParams & params )
 {
  m_correctionName = correctionName ;
  m_statusBits = statusBits ;
  m_chiSquare = chiSquare ;
  m_params = params ;
 }

void CalCorToolResult::Clear( Option_t * )
 {
  m_correctionName = "" ;
  m_statusBits = 0 ;
  m_chiSquare = 0. ;
  m_params.Clear() ;
 }

// dummy data, just for tests
void CalCorToolResult::Fake( UInt_t /* rank */, Float_t /* randNum */ )
 {
 }
 
// for tests
Bool_t CalCorToolResult::Compare( const CalCorToolResult & cl ) const {

    Bool_t result = true ;
    
//    result = result && rootdatautil::Compare(getEnergy(),cl.getEnergy(),"Energy") ;
//    result = result && rootdatautil::Compare(getPosition(),cl.getPosition(),"Position") ;
//    result = result && rootdatautil::Compare(getRmsSpread(),cl.getRmsSpread(),"RmsSpread") ;
//
//    if (!result) {
//        std::cout<<"Comparison ERROR for CalCorToolResult"<<std::endl ;
//    }
    return result ;

}

CalCorToolResult::CalCorToolResult()
 { Clear() ; }

CalCorToolResult::CalCorToolResult
  ( const TString & correctionName,
    UInt_t statusBits,
    Double_t chiSquare,
    const CalParams & params )
 { init(correctionName,statusBits,chiSquare,params) ; }

void CalCorToolResult::Print( Option_t * ) const
 {
//  std::cout
//    << m_energy
//    << " " << m_position.X() << "/" << m_position.Y() << "/" << m_position.Z()
//    << " " << m_rmsSpread.X() << "/" << m_rmsSpread.Y() << "/" << m_rmsSpread.Z() ;
 }


CalCorToolResult::~CalCorToolResult()
 {}

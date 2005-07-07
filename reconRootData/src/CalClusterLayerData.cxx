
#include <reconRootData/CalClusterLayerData.h>
#include <commonRootData/RootDataUtil.h>
#include <iostream>

ClassImp(CalClusterLayerData)

void CalClusterLayerData::init( Double_t energy, TVector3 position, TVector3 rmsSpread )
 {
  m_energy = energy ;
  m_position = position ;
  m_rmsSpread = rmsSpread ;
 }

void CalClusterLayerData::Clear( Option_t * )
 { init(0.,TVector3(0.,0.,0.),TVector3(0.,0.,0.)) ; }

// dummy data, just for tests
void CalClusterLayerData::Fake( UInt_t /* rank */, Float_t /* randNum */ )
 { init(1.,TVector3(2.,3.,4.),TVector3(5.,6.,7.)) ; }
 
// for tests
Bool_t CalClusterLayerData::Compare( const CalClusterLayerData & cl ) const {

    Bool_t result = true ;
    
    result = result && rootdatautil::Compare(getEnergy(),cl.getEnergy(),"Energy") ;
    result = result && rootdatautil::Compare(getPosition(),cl.getPosition(),"Position") ;
    result = result && rootdatautil::Compare(getRmsSpread(),cl.getRmsSpread(),"RmsSpread") ;

    if (!result) {
        std::cout<<"Comparison ERROR for CalClusterLayerData"<<std::endl ;
    }
    return result ;

}

CalClusterLayerData::CalClusterLayerData()
 { Clear() ; }

CalClusterLayerData::CalClusterLayerData( Double_t energy, TVector3 position, TVector3 rmsSpread )
 { init(energy,position,rmsSpread) ; }

void CalClusterLayerData::Print( Option_t * ) const
 {
  std::cout
    << m_energy
    << " " << m_position.X() << "/" << m_position.Y() << "/" << m_position.Z()
    << " " << m_rmsSpread.X() << "/" << m_rmsSpread.Y() << "/" << m_rmsSpread.Z() ;
 }

CalClusterLayerData::~CalClusterLayerData()
 {}

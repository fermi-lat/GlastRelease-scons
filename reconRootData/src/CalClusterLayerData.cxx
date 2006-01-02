
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
void CalClusterLayerData::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ )
 { init(1.,TVector3(2.,3.,4.),TVector3(5.,6.,7.)) ; }
 
// for tests
Bool_t CalClusterLayerData::CompareInRange( const CalClusterLayerData & cl, const std::string & name ) const {

    Bool_t result = true ;
    
    result = rootdatautil::CompareInRange(getEnergy(),cl.getEnergy(),"Energy") && result ;
    result = rootdatautil::CompareInRange(getPosition(),cl.getPosition(),"Position") && result ;
    result = rootdatautil::CompareInRange(getRmsSpread(),cl.getRmsSpread(),"RmsSpread") && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for CalClusterLayerData"<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
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

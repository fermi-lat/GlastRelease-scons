
#include <reconRootData/CalClusterLayerData.h>
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

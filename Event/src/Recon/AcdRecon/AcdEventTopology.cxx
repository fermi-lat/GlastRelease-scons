// File and Version information:
// $Header$
//
//  Implementation file of AcdEventTopology class
//  
// Authors:
//
//    Eric Charles
//
//

#include "Event/Recon/AcdRecon/AcdEventTopology.h"

#include "GaudiKernel/MsgStream.h"

using namespace Event;

AcdEventTopology::AcdEventTopology()
  :m_tileCount(0),
   m_ribbonCount(0),
   m_tileVeto(0),
   m_tileEnergy(0.),
   m_ribbonEnergy(0.),
   m_nTilesTop(0),
   m_nVetoTop(0),
   m_tileEnergyTop(0),
   m_nSidesHit(0),
   m_nSidesVeto(0){
  for ( unsigned i(0); i < 4; i++ ) {
    m_nTilesSideRow[i] = 0;
    m_nTilesSideFace[i] = 0;
    m_nVetoSideRow[i] = 0;
    m_nVetoSideFace[i] = 0;
    m_tileEnergySideRow[i] = 0.;
    m_tileEnergySideFace[i] = 0.;
  }
}

/// Copy constructor
AcdEventTopology::AcdEventTopology(const AcdEventTopology& other)
  :m_tileCount(other.m_tileCount),
   m_ribbonCount(other.m_ribbonCount),
   m_tileVeto(other.m_tileVeto),
   m_tileEnergy(other.m_tileEnergy),
   m_ribbonEnergy(other.m_ribbonEnergy),   
   m_nTilesTop(other.m_nTilesTop),
   m_nVetoTop(other.m_nVetoTop),
   m_tileEnergyTop(other.m_tileEnergyTop),
   m_nSidesHit(other.m_nSidesHit),
   m_nSidesVeto(other.m_nSidesVeto){
  for ( unsigned i(0); i < 4; i++ ) {
    m_nTilesSideRow[i] = other.m_nTilesSideRow[i];
    m_nTilesSideFace[i] = other.m_nTilesSideFace[i] ;
    m_nVetoSideRow[i] = other.m_nVetoSideRow[i];
    m_nVetoSideFace[i] = other.m_nVetoSideFace[i] ;
    m_tileEnergySideRow[i] = other.m_tileEnergySideRow[i];
    m_tileEnergySideFace[i] = other.m_tileEnergySideFace[i];
  }  
}


/// Constructor for use in reconstruction, 
AcdEventTopology::AcdEventTopology(unsigned tileCount, unsigned ribbonCount, unsigned tileVeto,
                                   float tileEnergy, float ribbonEnergy,
                                   unsigned nTilesTop, unsigned nTilesSideRow[4], unsigned nTilesSideFace[4],
                                   unsigned nVetoTop, unsigned nVetoSideRow[4], unsigned nVetoSideFace[4],
                                   float tileEnergyTop, float tileEnergySideRow[4], float tileEnergySideFace[4],
                                   unsigned nSidesHit, unsigned nSidesVeto) 
  :m_tileCount(tileCount),
   m_ribbonCount(ribbonCount),
   m_tileVeto(tileVeto),
   m_tileEnergy(tileEnergy),
   m_ribbonEnergy(ribbonEnergy),
   m_nTilesTop(nTilesTop),
   m_nVetoTop(nVetoTop),
   m_tileEnergyTop(tileEnergyTop),
   m_nSidesHit(nSidesHit),
   m_nSidesVeto(nSidesVeto){
  for ( unsigned i(0); i < 4; i++ ) {
    m_nTilesSideRow[i] = nTilesSideRow[i];
    m_nTilesSideFace[i] = nTilesSideFace[i] ;
    m_nVetoSideRow[i] = nVetoSideRow[i];
    m_nVetoSideFace[i] = nVetoSideFace[i] ;
    m_tileEnergySideRow[i] = tileEnergySideRow[i];
    m_tileEnergySideFace[i] = tileEnergySideFace[i];
  }
}

const AcdEventTopology& AcdEventTopology::operator=(const AcdEventTopology& other) {
  if ( this == &other ) { return *this; }
  m_tileCount = other.m_tileCount;
  m_ribbonCount = other.m_ribbonCount;
  m_tileVeto = other.m_tileVeto;
  m_tileEnergy = other.m_tileEnergy;
  m_ribbonEnergy = other.m_ribbonEnergy;
  m_nTilesTop = other.m_nTilesTop;
  m_tileEnergyTop = other.m_tileEnergyTop;
  m_nSidesHit = other.m_nSidesHit;
  m_nSidesVeto = other.m_nSidesVeto;
  for ( unsigned i(0); i < 4; i++ ) {
    m_nTilesSideRow[i] = other.m_nTilesSideRow[i];
    m_nTilesSideFace[i] = other.m_nTilesSideFace[i] ;
    m_nVetoSideRow[i] = other.m_nVetoSideRow[i];
    m_nVetoSideFace[i] = other.m_nVetoSideFace[i] ;
    m_tileEnergySideRow[i] = other.m_tileEnergySideRow[i];
    m_tileEnergySideFace[i] = other.m_tileEnergySideFace[i];
  }  
  return *this;
}

void AcdEventTopology::set(unsigned tileCount, unsigned ribbonCount, unsigned tileVeto,
                           float tileEnergy, float ribbonEnergy,
                           unsigned nTilesTop, unsigned nTilesSideRow[4], unsigned nTilesSideFace[4],
                           unsigned nVetoTop, unsigned nVetoSideRow[4], unsigned nVetoSideFace[4],
                           float tileEnergyTop, float tileEnergySideRow[4], float tileEnergySideFace[4],
                           unsigned nSidesHit, unsigned nSidesVeto) {
  m_tileCount = tileCount;
  m_ribbonCount = ribbonCount;
  m_tileVeto = tileVeto;
  m_tileEnergy = tileEnergy;
  m_ribbonEnergy = ribbonEnergy;
  m_nTilesTop = nTilesTop;
  m_nVetoTop = nVetoTop;
  m_tileEnergyTop = tileEnergyTop;
  m_nSidesHit = nSidesHit;
  m_nSidesVeto = nSidesVeto;
  for ( unsigned i(0); i < 4; i++ ) {
    m_nTilesSideRow[i] = nTilesSideRow[i];
    m_nTilesSideFace[i] = nTilesSideFace[i] ;
    m_nVetoSideRow[i] = nVetoSideRow[i];
    m_nVetoSideFace[i] = nVetoSideFace[i] ;
    m_tileEnergySideRow[i] = tileEnergySideRow[i];
    m_tileEnergySideFace[i] = tileEnergySideFace[i];
  }
}

void AcdEventTopology::writeOut(MsgStream& stream) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << MSG::DEBUG
         << "AcdEventTopology. Tiles: " << m_tileCount << ", Ribbons: " << m_ribbonCount << " ,Vetos: " << m_tileVeto << std::endl
         << ", E_tile: " << m_tileEnergy << ", E_rib: " << m_ribbonEnergy << std::endl
         << "NTtile.  Top: " << m_nTilesTop 
         << ", SideRows: " << m_nTilesSideRow[0] << ' ' << m_nTilesSideRow[1] << ' ' << m_nTilesSideRow[2] << ' ' << m_nTilesSideRow[3]
         << ", Faces: " << m_nTilesSideFace[0] << ' ' << m_nTilesSideFace[1] << ' ' << m_nTilesSideFace[2] << ' ' << m_nTilesSideFace[3] << std::endl
         << "NVeto.  Top: " << m_nVetoTop 
         << ", SideRows: " << m_nVetoSideRow[0] << ' ' << m_nVetoSideRow[1] << ' ' << m_nVetoSideRow[2] << ' ' << m_nVetoSideRow[3]
         << ", Faces: " << m_nVetoSideFace[0] << ' ' << m_nVetoSideFace[1] << ' ' << m_nVetoSideFace[2] << ' ' << m_nVetoSideFace[3] << std::endl
         << "Energy.  Top: " << m_tileEnergyTop
         << ", SideRows: " << m_tileEnergySideRow[0] << ' ' << m_tileEnergySideRow[1] << ' ' 
         << m_tileEnergySideRow[2] << ' ' << m_tileEnergySideRow[3] 
         << ", Faces: " << m_tileEnergySideFace[0] << ' ' << m_tileEnergySideFace[1] << ' ' 
         << m_tileEnergySideFace[2] << ' ' << m_tileEnergySideFace[3] << std::endl
         << "SidesHit: " << m_nSidesHit << ", SideVetoed: " << m_nSidesVeto 
         << endreq;
}



void AcdEventTopology::ini()
// Purpose: reset all data members to 0
//
{ 
  m_tileCount = 0;
  m_ribbonCount = 0;
  m_tileVeto = 0;
  m_tileEnergy = 0.;
  m_ribbonEnergy = 0.;
  m_nVetoTop = 0;
  m_nTilesTop = 0;
  m_tileEnergyTop = 0;
  m_nSidesHit = 0;
  m_nSidesVeto = 0;
  for ( unsigned i(0); i < 4; i++ ) {
    m_nTilesSideRow[i] = 0;
    m_nTilesSideFace[i] = 0;
    m_nVetoSideRow[i] = 0;
    m_nVetoSideFace[i] = 0;
    m_tileEnergySideRow[i] = 0.;
    m_tileEnergySideFace[i] = 0.;
  }
}
 

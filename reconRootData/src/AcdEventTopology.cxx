#include "reconRootData/AcdEventTopology.h"
#include "Riostream.h"

ClassImp(AcdEventTopology) ;

/// Default constructor.  Set everything to default values
AcdEventTopology::AcdEventTopology() 
  :TObject(){
  Clear("");
}

/// Copy constructor
AcdEventTopology::AcdEventTopology(const AcdEventTopology& other)
  :TObject(other),
   m_tileCount(other.m_tileCount),
   m_ribbonCount(other.m_ribbonCount),
   m_vetoCount(other.m_vetoCount),
   m_tileVeto(other.m_tileVeto),
   m_totalTileEnergy(other.m_totalTileEnergy),
   m_totalRibbonEnergy(other.m_totalRibbonEnergy),
   m_tileEnergy(other.m_tileEnergy),
   m_ribbonEnergy(other.m_ribbonEnergy),
   m_ghostTileEnergy(other.m_ghostTileEnergy),
   m_ghostRibbonEnergy(other.m_ghostRibbonEnergy),
   m_triggerTileEnergy(other.m_triggerTileEnergy),
   m_triggerRibbonEnergy(other.m_triggerRibbonEnergy),
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

/// Assignment operator
const AcdEventTopology& AcdEventTopology::operator=(const AcdEventTopology& other)
{
  if ( this == &other ) return *this;
  m_tileCount = other.m_tileCount;
  m_ribbonCount = other.m_ribbonCount;
  m_vetoCount = other.m_vetoCount;
  m_tileVeto = other.m_tileVeto;
  m_totalTileEnergy = other.m_totalTileEnergy;
  m_totalRibbonEnergy = other.m_totalRibbonEnergy;
  m_tileEnergy = other.m_tileEnergy;
  m_ribbonEnergy = other.m_ribbonEnergy;
  m_ghostTileEnergy = other.m_ghostTileEnergy;
  m_ghostRibbonEnergy = other.m_ghostRibbonEnergy;
  m_triggerTileEnergy = other.m_triggerTileEnergy;
  m_triggerRibbonEnergy = other.m_triggerRibbonEnergy;
  m_nTilesTop = other.m_nTilesTop;
  m_nVetoTop = other.m_nVetoTop;
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

/// set all the values
void AcdEventTopology::set(UInt_t tileCount, UInt_t ribbonCount, UInt_t vetoCount, UInt_t tileVeto,                                    
                           Float_t totalTileEnergy, Float_t totalRibbonEnergy, Float_t tileEnergy, Float_t ribbonEnergy,               
                           Float_t ghostTileEnergy, Float_t ghostRibbonEnergy, Float_t triggerTileEnergy, Float_t triggerRibbonEnergy, 
                           UInt_t nTilesTop, UInt_t nTilesSideRow[4], UInt_t nTilesSideFace[4],                                        
                           UInt_t nVetoTop, UInt_t nVetoSideRow[4], UInt_t nVetoSideFace[4],                                           
                           Float_t tileEnergyTop, Float_t tileEnergySideRow[4], Float_t tileEnergySideFace[4],                         
                           UInt_t nSidesHit, UInt_t nSidesVeto) {
  m_tileCount = tileCount;
  m_ribbonCount = ribbonCount;
  m_vetoCount = vetoCount;
  m_tileVeto = tileVeto;
  m_totalTileEnergy = totalTileEnergy;
  m_totalRibbonEnergy = totalRibbonEnergy;
  m_tileEnergy = tileEnergy;
  m_ribbonEnergy = ribbonEnergy;
  m_ghostTileEnergy = ghostTileEnergy;
  m_ghostRibbonEnergy = ghostRibbonEnergy;
  m_triggerTileEnergy = triggerTileEnergy;
  m_triggerRibbonEnergy = triggerRibbonEnergy;
  m_nTilesTop = nTilesTop;
  m_nVetoTop = nVetoTop;
  m_tileEnergyTop = tileEnergyTop;
  m_nSidesHit = nSidesHit;
  m_nSidesVeto = nSidesVeto;
  for ( UInt_t i(0); i < 4; i++ ) {
    m_nTilesSideRow[i] = nTilesSideRow[i];
    m_nTilesSideFace[i] = nTilesSideFace[i] ;
    m_nVetoSideRow[i] = nVetoSideRow[i];
    m_nVetoSideFace[i] = nVetoSideFace[i] ;
    m_tileEnergySideRow[i] = tileEnergySideRow[i];
    m_tileEnergySideFace[i] = tileEnergySideFace[i];
  }
}  

/// reset all the values to their default
void AcdEventTopology::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_tileCount = 0;
  m_ribbonCount = 0;
  m_vetoCount = 0;
  m_tileVeto = 0;
  m_totalTileEnergy = 0.;
  m_totalRibbonEnergy = 0.;
  m_tileEnergy = 0.;
  m_ribbonEnergy = 0.;
  m_ghostTileEnergy = 0.;
  m_ghostRibbonEnergy = 0.;
  m_triggerTileEnergy = 0.;
  m_triggerRibbonEnergy = 0.;
  m_nTilesTop = 0;
  m_nVetoTop = 0;
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

/// Print out this structure
void AcdEventTopology::Print(Option_t *option) const {
  std::cout  << "AcdEventTopology." 
             << " Tiles: " << m_tileCount << ", Ribbons: " << m_ribbonCount << " ,Vetos: " << m_tileVeto << std::endl
             << ", E_tile: " << m_totalTileEnergy << ", E_rib: " << m_totalRibbonEnergy << std::endl
             << "NTtile.  Top: " << m_nTilesTop 
             << ", SideRows: " << m_nTilesSideRow[0] << ' ' << m_nTilesSideRow[1] << ' ' << m_nTilesSideRow[2] << ' ' << m_nTilesSideRow[3]
             << ", Faces: " << m_nTilesSideFace[0] << ' ' << m_nTilesSideFace[1] << ' ' << m_nTilesSideFace[2] << ' ' << m_nTilesSideFace[3] << std::endl
             << "NTtile.  Top: " << m_nVetoTop 
             << ", SideRows: " << m_nVetoSideRow[0] << ' ' << m_nVetoSideRow[1] << ' ' << m_nVetoSideRow[2] << ' ' << m_nVetoSideRow[3]
             << ", Faces: " << m_nVetoSideFace[0] << ' ' << m_nVetoSideFace[1] << ' ' << m_nVetoSideFace[2] << ' ' << m_nVetoSideFace[3] << std::endl
             << "Energy.  Top: " << m_tileEnergyTop
             << ", SideRows: " << m_tileEnergySideRow[0] << ' ' << m_tileEnergySideRow[1] << ' ' 
             << m_tileEnergySideRow[2] << ' ' << m_tileEnergySideRow[3] 
             << ", Faces: " << m_tileEnergySideFace[0] << ' ' << m_tileEnergySideFace[1] << ' ' 
             << m_tileEnergySideFace[2] << ' ' << m_tileEnergySideFace[3] << std::endl
             << "SidesHit: " << m_nSidesHit << ", SideVetoed: " << m_nSidesVeto 
             << std::endl;
}

void AcdEventTopology::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdEventTopology::CompareInRange( const AcdEventTopology& /* other */, const std::string & /* name */ ) const {
  return kTRUE;
}

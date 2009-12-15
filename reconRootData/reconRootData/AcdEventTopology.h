#ifndef ROOT_AcdEventTopology_H
#define ROOT_AcdEventTopology_H

#include "TObject.h"
#include <string>

/** 
 * @class AcdEventTopology
 * @brief ROOT object which stores information about the Point of Closest Approach (POCA) between an 
 * extrapolated track and a hit Acd element (tile or ribbon).  
 *  
 * Most of the structure of the object comes from the base classes AcdTkrLocalCoords and AcdPocaData
 * 
 * The class adds only enough information to define the involved in the POCA.
 *    - const idents::AcdId& getId()  
 *      - which returns the ID of the hit element
 *    - const int getTrackIndex()  
 *      - which returns the index of the track which did the hitting
 *
 * @author Eric Charles
 *
 * $Header$
 */



class AcdEventTopology : public TObject {

public:

  /// Default constructor.  
  AcdEventTopology();
  
  /// Copy constructor
  AcdEventTopology(const AcdEventTopology& params);

  /// Assignment operator
  const AcdEventTopology& operator=(const AcdEventTopology& params);

  /// D'tor
  virtual ~AcdEventTopology(){}

    /// Direct access to parameters
  inline UInt_t getTileCount() const { return m_tileCount; }
    
  inline UInt_t getRibbonCount() const { return m_ribbonCount; }

  inline UInt_t getTileVeto() const { return m_tileVeto; }

  inline Float_t getTileEnergy() const { return m_tileEnergy; }

  inline Float_t getRibbonEnergy() const { return m_ribbonEnergy; }

  inline UInt_t getVetoCountTop() const { return m_nVetoTop; }
  
  inline UInt_t getVetoCountSideRow(unsigned iRow) const { return m_nVetoSideRow[iRow]; }
  
  inline UInt_t getVetoCountSideFace(unsigned iFace) const { return m_nVetoSideFace[iFace-1]; }

  inline UInt_t getTileCountTop() const { return m_nTilesTop; }

  inline UInt_t getTileCountSideRow(UInt_t iRow) const { return m_nTilesSideRow[iRow]; }

  inline UInt_t getTileCountSideFace(UInt_t iFace) const { return m_nTilesSideFace[iFace-1]; }

  inline Float_t getTileEnergyTop() const { return m_tileEnergyTop; }

  inline Float_t getTileEnergySideRow(UInt_t iRow) const { return m_tileEnergySideRow[iRow]; }

  inline Float_t getTileEnergySideFace(UInt_t iFace) const { return m_tileEnergySideFace[iFace-1]; }

  inline UInt_t getNSidesHit() const { return m_nSidesHit; }

  inline UInt_t getNSidesVeto() const { return m_nSidesVeto; }

  /// set all the values
  void set(UInt_t tileCount, UInt_t ribbonCount, UInt_t tileVeto,
	   Float_t tileEnergy, Float_t ribbonEnergy,
	   UInt_t nTilesTop, UInt_t nTilesSideRow[4], UInt_t nTilesSideFace[4],
	   UInt_t nVetoTop, UInt_t nVetoSideRow[4], UInt_t nVetoSideFace[4],
	   Float_t tilesEnergyTop, Float_t tileEnergySideRow[4], Float_t tileEnergySideFace[4],
	   UInt_t nSidesHit, UInt_t nSidesVeto);
  
  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdEventTopology&, const std::string& name = "" ) const ; // for tests   

private:
    
  UInt_t m_tileCount;
  
  UInt_t m_ribbonCount;
  
  UInt_t m_tileVeto;
  
  Float_t m_tileEnergy;
  
  Float_t m_ribbonEnergy;
  
  UInt_t m_nTilesTop;
  
  UInt_t m_nTilesSideRow[4];
  
  UInt_t m_nTilesSideFace[4];
  
  UInt_t m_nVetoTop;
  
  UInt_t m_nVetoSideRow[4];
  
  UInt_t m_nVetoSideFace[4];

  Float_t m_tileEnergyTop;
  
  Float_t m_tileEnergySideRow[4];
  
  Float_t m_tileEnergySideFace[4];
  
  UInt_t m_nSidesHit;

  UInt_t m_nSidesVeto;
  
  ClassDef(AcdEventTopology,1)
  
};
#endif

// $Header$

#ifndef INTEGRATINGDETECTORMANAGER_H
#define INTEGRATINGDETECTORMANAGER_H


#include "DisplayManager.h"

#include "idents/VolumeIdentifier.h"

#include "GlastEvent/MonteCarlo/McIntegratingHit.h"

#include "DetectorManager.h"
#include <map>

class DetectorConstruction;
class IDataProviderSvc;

/** 
 * @class IntDetectorManager
 *
 * @brief A concrete DetectorManager
 *
 * This class implement the abstract DetectorManager; this is used for all the
 * detector sensitive volume that need to register hits as McIntegratingHits.
 * 
 * @author T.Burnett and R.Giannitrapani
 *    
 * $Header$
 */
class IntDetectorManager : public DetectorManager {
public:
    
  //! @param det the DetectorConstruction pointer to retrive the map of volume
  //!        ids for all sensitive detectors 
  //! @param esv the data provider service for TDS access 
  IntDetectorManager( DetectorConstruction*, IDataProviderSvc*);
  
  //! Clears things; this implement a pure abstract method in the
  //! hierarchy ancestor of this class (geant4 name convention)
  virtual void Initialize(G4HCofThisEvent*);
  
  //! Called by G4 in each step in a sensitive volume; this implement a pure
  //! abstract method in the hierarchy ancestor of this class (geant4 name
  //! convention)
  virtual G4bool ProcessHits(G4Step* aStep ,G4TouchableHistory*);
  
  //! End of event will finish hits retrival; this implement a pure abstract
  //! method in the hierarchy ancestor of this class (geant4 name convention)
  virtual void EndOfEvent(G4HCofThisEvent*);
  
 private:
  /// The collection of McIntegratingHit to save in the TDS
  McIntegratingHitVector *m_intHit;  
  
  /// A map of McIntegratingHit indicized by volume id to pile up energy
  /// deposited
  std::map<idents::VolumeIdentifier,mc::McIntegratingHit*> m_detectorList;
};
#endif

#ifndef POSITIONDETECTORMANAGER_H
#define POSITIONDETECTORMANAGER_H

#include "Event/MonteCarlo/McPositionHit.h"

#include "DetectorManager.h"

class DetectorConstruction;
class IDataProviderSvc;

/** 
 * @class PosDetectorManager
 *
 * @brief A concrete DetectorManager
 *
 * This class implement the abstract DetectorManager; this is used for all the
 * detector sensitive volume that need to register hits as McPositionHits.
 * 
 * @author T.Burnett and R.Giannitrapani
 *    
 * $Header$
 */
class PosDetectorManager : public DetectorManager {
public:

  /// @param det the DetectorConstruction pointer to retrive the map of volume
  ///        ids for all sensitive detectors 
  /// @param esv the data provider service for TDS access 
  PosDetectorManager( DetectorConstruction* det, IDataProviderSvc* esv);
    
  /// Clears things; this implement a pure abstract method in the
  /// hierarchy ancestor of this class (geant4 name convention)
  virtual void Initialize(G4HCofThisEvent*);
  
  /// Called by G4 in each step in a sensitive volume; this implement a pure
  /// abstract method in the hierarchy ancestor of this class (geant4 name
  /// convention)
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  
  /// End of event will finish hits retrival; this implement a pure abstract
  /// method in the hierarchy ancestor of this class (geant4 name convention)
  virtual void EndOfEvent(G4HCofThisEvent*);
  
 private:
  /// The collection of McPositionHit to save in the TDS
  Event::McPositionHitVector *m_posHit;  
};
#endif

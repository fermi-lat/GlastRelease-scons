// $Header$

#ifndef POSITIONDETECTORMANAGER_H
#define POSITIONDETECTORMANAGER_H




#include "GlastEvent/MonteCarlo/McPositionHit.h"

#include "DetectorManager.h"

class DetectorConstruction;
class IDataProviderSvc;

/**

  */
class PosDetectorManager : public DetectorManager {
public:
    
    //! constructor called with pointer to DetectorConstruction, for map of (partial) ids
    //! needed for constructing the from the list of physical volumes in the touchable history
    //! @param idmap map of volume ids for all sensitive detectors
    PosDetectorManager( DetectorConstruction*, IDataProviderSvc*);
    
    //! initialize clears things 
    virtual void Initialize(G4HCofThisEvent*);

    //! G4 passes in each step in a sensitive volume
    virtual G4bool ProcessHits(G4Step* aStep ,G4TouchableHistory*);

    //! End of event will finish digitization
    virtual void EndOfEvent(G4HCofThisEvent*);
    
private:
    /// The collection of McPositionHit to save in the TDS
    McPositionHitVector *m_posHit;  

};

#endif

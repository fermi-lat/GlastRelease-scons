// $Header$

#ifndef POSITIONDETECTORMANAGER_H
#define POSITIONDETECTORMANAGER_H
#ifdef WIN32 // for G4 
#include <float.h>
#endif

#include "G4LogicalVolume.hh"
#include "G4VSensitiveDetector.hh"

#include "DisplayManager.h"

#include "idents/VolumeIdentifier.h"
#include "DetectorConstruction.h"

#include "GlastEvent/MonteCarlo/McPositionHit.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "DetectorManager.h"
#include <map>

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

    typedef std::map<idents::VolumeIdentifier, unsigned int> DetectorList;
    DetectorList m_detectorList;
};

#endif

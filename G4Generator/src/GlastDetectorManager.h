// $Header$

#ifndef GlastDetectorManager_h
#define GlastDetectorManager_h
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

#include <map>

/**
A single class that manages the GlastDetector hierarchy in the G4 environment:
- Reads the xml description file and instantiates them
- connects each GlastDetector object with a special interface object to convert G4Step information
- When called with a G4Step, construcs the ID from the touchable history, then looks up the 
  corresponding object

  */
class GlastDetectorManager : public G4VSensitiveDetector {
public:
    
    //! constructor called with pointer to DetectorConstruction, for map of (partial) ids
    //! needed for constructing the from the list of physical volumes in the touchable history
    //! @param idmap map of volume ids for all sensitive detectors
    GlastDetectorManager( DetectorConstruction*, IDataProviderSvc*);
    
    ~GlastDetectorManager();
    
    //! initialize clears things 
    virtual void Initialize(G4HCofThisEvent*);

    //! G4 passes in each step in a sensitive volume
    virtual G4bool ProcessHits(G4Step* aStep ,G4TouchableHistory*);

    //! End of event will finish digitization
    virtual void EndOfEvent(G4HCofThisEvent*);
    
    idents::VolumeIdentifier constructId(G4Step * aStep);
    
    //! Called from DetectorConstruction to set the sensitive detector propery
    void process(G4LogicalVolume*);

private:
    DetectorConstruction::IdMap* m_idMap;
    
    //! keep track of detectors
    typedef std::map<idents::VolumeIdentifier, unsigned int> DetectorList;
    typedef std::map<idents::VolumeIdentifier, double>DetectorEnergyTotal;
    DetectorList m_detectorList;
    DetectorEnergyTotal m_detectorEnergy;

    /// The pointer to the IdataProviderSvc
    IDataProviderSvc* m_esv;
    /// The collection of McPositionHit to save in the TDS
    McPositionHitVector *m_posHit;  
    
};

#endif
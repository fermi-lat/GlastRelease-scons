// $Header$

#ifndef DetectorManager_h
#define DetectorManager_h
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
class G4TouchableHistory;
namespace mc {class McPositionHit;}

/**
A single class that manages the GlastDetector hierarchy in the G4 environment:
- Reads the xml description file and instantiates them
- connects each GlastDetector object with a special interface object to convert G4Step information
- When called with a G4Step, construcs the ID from the touchable history, then looks up the 
  corresponding object

  */
class DetectorManager : public G4VSensitiveDetector {
public:
    
    //! constructor called with pointer to DetectorConstruction, for map of (partial) ids
    //! needed for constructing the from the list of physical volumes in the touchable history
    //! @param idmap map of volume ids for all sensitive detectors
    DetectorManager( DetectorConstruction*, IDataProviderSvc* , std::string name);
    
    ~DetectorManager();
    //! Called from DetectorConstruction to set the sensitive detector propery
    void process(G4LogicalVolume*);
protected:

    idents::VolumeIdentifier constructId(G4Step * aStep);

    void makeBox(G4TouchableHistory* touched);

    /// The pointer to the IdataProviderSvc
    IDataProviderSvc* m_esv;
    
    void display(G4TouchableHistory* touched, mc::McPositionHit * hit);

private:
    DetectorConstruction::IdMap* m_idMap;

    //! keep track of hit detectors for display
    typedef std::map<idents::VolumeIdentifier, unsigned int> DetectorList;
    DetectorList m_detectorList;

#if 0
    
    typedef std::map<idents::VolumeIdentifier, double>DetectorEnergyTotal;
    DetectorEnergyTotal m_detectorEnergy;
    /// The collection of McPositionHit to save in the TDS
    McPositionHitVector *m_posHit;  
    
#endif
};

#endif

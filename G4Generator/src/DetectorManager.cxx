// File and Version Information:
// $Header$
//
// Description: This is an abstract class that represent a generic sensitive
// detectors manager. Its methods are common to both kind of detectors managers
// of Glast, i.e. PosDetectorManager and IntDetectorManager
//
// Author(s):
//      R.Giannitrapani

#include "DetectorManager.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "idents/VolumeIdentifier.h"

// Geant4 interface
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"

DetectorManager::DetectorManager(DetectorConstruction::IdMap *map,
                                 IDataProviderSvc* esv, IGlastDetSvc* gsv, std::string name)
  :m_idMap(map), m_esv(esv), m_gsv(gsv), G4VSensitiveDetector(name)
{
  // Inputs: the IdMap, coming from the DetectorConstruction and actually built
  // from the GlastDetSvc calls to detModel functionalities, provide a mapping
  // between Geant4 physical volumes and VolumeIdentifiers. The DataProviderSvc
  // is used by subclasses to save hits objects in the TDS. The name is used to
  // distinguish different detectors managers by the SDManager of Geant4.
  // Restrictions and Caveats: Note that this is an abstract class, so it cannot
  // be instanciated by itself
  
  // tell G4 about us
  G4SDManager::GetSDMpointer()->AddNewDetector( this );
}

DetectorManager::~DetectorManager()
{
}
   
idents::VolumeIdentifier DetectorManager::constructId(G4Step * aStep)
{
  // Purpouse and Methods: this method build a VolumeIdentifier by using the
  // IdMap and information in the aStep; in practice it concatenates identifiers
  // associated, with the help of the IdMap, to the hierarchy of touchable
  // volumes starting from the volume containing the aStep and ending in the
  // world volume
  // Inputs: the G4Step pointer aStep
  using  idents::VolumeIdentifier;
  VolumeIdentifier ret;
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) {
    const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
//**    if( physVol->GetMother()==0) break;
    VolumeIdentifier id = (*m_idMap)[physVol];
    ret.prepend(id);
  }
  return ret;
}

void DetectorManager::process(G4LogicalVolume* lvol)
{
  // Purpouse and Mehtods: this method is used to register a volume (flagged in
  // the xml description as sensitive) to this manager
  // Inputs: the logical volume pointer
  lvol->SetSensitiveDetector(this);
}


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
                                 IDataProviderSvc* esv, std::string name)
  :m_idMap(map), m_esv(esv),G4VSensitiveDetector(name)
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
  // this will be null if no GUI
  m_display = DisplayManager::instance();
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
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) {
    const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
    if( physVol->GetMother()==0) break;
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


void DetectorManager::makeDisplayBox(G4TouchableHistory* touched,
                                     idents::VolumeIdentifier id,
                                     bool hitBox)
{
  // Purpouse and Methods: this method draw a box in the gui (if it has been
  // activated) by recovering the global 3d transformation from the volume
  // touchables hierarchy 
  // Inputs: the hierarchy of touchable volumes touched, the identifier id (to
  // be shown along the box in the GUI) and a flag that says if the box is a
  // position one or an integrating one

  if(displayMgr()==0)return;
  G4VPhysicalVolume* pvol = touched->GetVolume(); 
    
  HepTransform3D 
    global(*(touched->GetRotation()), 
           touched->GetTranslation());
    
    
  const G4LogicalVolume* lvol = pvol->GetLogicalVolume();
  const G4VSolid * solid = lvol->GetSolid();
  const G4Box* box = dynamic_cast<const G4Box*>(solid);
  if( box !=0){
    double 
      x = 2*box->GetXHalfLength(), 
      y = 2*box->GetYHalfLength(), 
      z = 2*box->GetZHalfLength();
        
    if( hitBox)  DisplayManager::instance()->addHitBox(global, x,y,z);
    else DisplayManager::instance()->addIntegratingBox(global, x,y,z);


  }
  displayMgr()->addIdDisplay(global, id);
    
}

void DetectorManager::display(G4TouchableHistory* touched, 
                              idents::VolumeIdentifier id, 
                              const HepPoint3D& entry, const HepPoint3D& exit)
{
  // Purpouse and Methods: this method is used to draw on the GUI (if activated)
  // both a box and the step hit
  // Inputs: the touchable history touched, the volume identifier id and the
  // entry and exit points of the step

  if( displayMgr()==0) return;
  displayMgr()->addHit(entry, exit);

  if( m_detectorList[id]==0) {
    makeDisplayBox( touched , id, true);        
  }
  ++ m_detectorList[id]; 

}

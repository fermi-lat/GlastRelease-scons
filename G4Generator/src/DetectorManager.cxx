// $Header$

#include "DetectorManager.h"
#include <iostream>
#include "CLHEP/Geometry/Transform3D.h"
#include "GlastEvent/MonteCarlo/McPositionHit.h"
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
    // and tell G4 about us
    G4SDManager::GetSDMpointer()->AddNewDetector( this );
    // this will be null if no GUI
    m_display = DisplayManager::instance();
}

DetectorManager::~DetectorManager()
{
}
   
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idents::VolumeIdentifier DetectorManager::constructId(G4Step * aStep)
{
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
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DetectorManager::process(G4LogicalVolume* lvol)
{
    lvol->SetSensitiveDetector(this);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DetectorManager::makeDisplayBox(G4TouchableHistory* touched,
                                     idents::VolumeIdentifier id,
                                     bool hitBox)
{
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
        else        DisplayManager::instance()->addIntegratingBox(global, x,y,z);


    }
    displayMgr()->addIdDisplay(global, id);
    
}

void DetectorManager::display(G4TouchableHistory* touched, 
                              idents::VolumeIdentifier id, 
                              const HepPoint3D& entry, const HepPoint3D& exit)
{
    if( displayMgr()==0) return;
    displayMgr()->addHit(entry, exit);

    if( m_detectorList[id]==0) {
        makeDisplayBox( touched , id, true);        
    }
    ++ m_detectorList[id]; 

}
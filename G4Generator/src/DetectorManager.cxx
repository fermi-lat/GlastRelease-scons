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

DetectorManager::DetectorManager(DetectorConstruction *det,
                                 IDataProviderSvc* esv, std::string name)
:m_esv(esv),G4VSensitiveDetector(name)
{
    // Recover the physicals-id map from the detectorconstruction
    m_idMap = det->idMap();

    // and tell G4 about us
    G4SDManager::GetSDMpointer()->AddNewDetector( this );

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
void DetectorManager::makeBox(G4TouchableHistory* touched)
{
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
        
        DisplayManager::instance()->addHitBox(global, x,y,z);
    }
    
}
void DetectorManager::display(G4TouchableHistory* touched, mc::McPositionHit * hit)
{
    
    DisplayManager::instance()->addHit(hit->entryPoint(), hit->exitPoint());

    if( m_detectorList[hit->volumeID()]==0) {
        makeBox( touched );        
    }
    ++ m_detectorList[hit->volumeID()]; 

}
// $Header$

#include "PosDetectorManager.h"
#include <iostream>
#include "CLHEP/Geometry/Transform3D.h"
#include "GlastEvent/MonteCarlo/McPositionHit.h"
#include "idents/VolumeIdentifier.h"

// Geant4 interface
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"


namespace{
    void makeBox(G4TouchableHistory* touched)
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
}

PosDetectorManager::PosDetectorManager(DetectorConstruction *det,
                                           IDataProviderSvc* esv)
:DetectorManager(det, esv,"PositionDetectorManager")
{
}

void PosDetectorManager::Initialize(G4HCofThisEvent*HCE)
{
  m_detectorList.clear();
  // At the start of the event we create a new container
  m_posHit = new McPositionHitVector;    
}

G4bool PosDetectorManager::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
    
    // Energy Deposition & Step Length
    
    G4double edep = aStep->GetTotalEnergyDeposit()/MeV;
    G4double stepl = aStep->GetStepLength()/mm;
    
    if ((edep==0.)) return false;          
    // Physical Volume
    
    G4TouchableHistory* theTouchable
        = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
    G4LogicalVolume* logVol = physVol->GetLogicalVolume();
    G4String material = logVol->GetMaterial()->GetName();
    G4String nameVolume = physVol->GetName();
    
    G4ThreeVector InitPos = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector FinPos = aStep->GetPostStepPoint()->GetPosition();
    
    // determine the ID by studying the history, then call appropriate 
    idents::VolumeIdentifier id = constructId(aStep);

    // Filling of the hits container
    mc::McPositionHit *hit = new mc::McPositionHit;
    
    hit->setDepositedEnergy(edep);
    hit->setVolumeID(id);
    hit->setEntryPoint(InitPos);
    hit->setExitPoint(FinPos);

    m_posHit->push_back(hit);

    //**** interface to display *************
    
    DisplayManager::instance()->addHit(InitPos, FinPos);
    
    if( m_detectorList[id]==0) {
        makeBox( theTouchable );        
    }
    ++ m_detectorList[id]; 

    return true;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PosDetectorManager::EndOfEvent(G4HCofThisEvent* HCE)
{
    // store the hits in the TDS
    m_esv->registerObject("/Event/MC/PositionHitsCol", m_posHit);    

    std::cout << "Actual Event done! " << m_posHit->size() 
        << " position hits stored in the TDS" << std::endl;

}


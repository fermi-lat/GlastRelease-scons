// $Header$

#include "IntDetectorManager.h"
#include <iostream>
#include "CLHEP/Geometry/Transform3D.h"
#include "GlastEvent/MonteCarlo/McIntegratingHit.h"
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

IntDetectorManager::IntDetectorManager(DetectorConstruction *det,
                                           IDataProviderSvc* esv)
:DetectorManager(det, esv,"IntegratingDetectorManager")
{
}

void IntDetectorManager::Initialize(G4HCofThisEvent*HCE)
{
  m_detectorList.clear();
  // At the start of the event we create a new container
  m_intHit = new McIntegratingHitVector;    
}

G4bool IntDetectorManager::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
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

    //**** interface to display *************
    
    //    DisplayManager::instance()->addHit(InitPos, FinPos);
    
    if( !m_detectorList[id]) {
        makeBox( theTouchable );        
    }

    
    // We fill an integrating hit
    mc::McIntegratingHit *hit; 
    // If the hit has already been created we use it, otherwise we
    // create a new one
    if( !(hit = m_detectorList[id]))
      {
	// Filling of the hits container
	hit = new mc::McIntegratingHit;
	m_intHit->push_back(hit);
	m_detectorList[id] = hit;
      }
	
    //    hit->setDepositedEnergy(edep);
    //	hit->setVolumeID(id);
    //	hit->setEntryPoint(InitPos);
    //	hit->setExitPoint(FinPos);



    return true;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void IntDetectorManager::EndOfEvent(G4HCofThisEvent* HCE)
{
    // store the hits in the TDS
    m_esv->registerObject("/Event/MC/IntegratingHitsCol", m_intHit);    

    std::cout << "Actual Event done! " << m_intHit->size() 
        << " integrating hits stored in the TDS" << std::endl;

}


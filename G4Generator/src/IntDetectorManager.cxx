// File and Version Information:
// $Header$
//
// Description: This is a concrete implementation of the DetectorManager
// abstract class; this one is used to manage sensitive detectors of integrating
// type
//
// Author(s):
//      R.Giannitrapani

#include "IntDetectorManager.h"
#include "McParticleManager.h"

#include <iostream>
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h"

// Geant4 interface
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"

#include <algorithm>

IntDetectorManager::IntDetectorManager(DetectorConstruction *det,
                                       IDataProviderSvc* esv)
  :DetectorManager(det->idMap(), esv,"IntegratingDetectorManager")
{
  // See the father class DetectorManager
}

void IntDetectorManager::Initialize(G4HCofThisEvent*)
{
  // Purpose and Method: initialization of the detectors manager
  // Inputs: the G4HCofThisEvent is hinerited from the Geant4 structure and is
  // of no use for our actual implementation of G4Generator (but must be there
  // for Geant4 internal working) 

  // clear the list of hited detectors
  m_detectorList.clear();
  // At the start of the event we create a new container for the TDS
  m_intHit = new Event::McIntegratingHitVector;    
}

G4bool IntDetectorManager::ProcessHits(G4Step* aStep,
                                       G4TouchableHistory* ROhist)
{    
  // Purpose and Method: this method is called internally by Geant4 every time a
  // particle issue an hit in a sensitive detector of this kind
  // Inputs: the step of the hit and the hierarchy of touchable volumes
  // Outpus: false if there is no hit to register, true otherwise

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

  // start position of the hit and final one
  G4ThreeVector prePos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector postPos = aStep->GetPostStepPoint()->GetPosition();
    
  // determine the ID by studying the history, then call appropriate 
  idents::VolumeIdentifier id = constructId(aStep);

  // We want to fill an integrating hit
  Event::McIntegratingHit *hit; 

  // If the hit has already been created we use it, otherwise we
  // create a new one
  if( !(hit = m_detectorList[id]))
    {
      // This draw the volume
      makeDisplayBox( theTouchable, id );        
      // A new object is needed
      hit = new Event::McIntegratingHit;
      // Set its volume identifier
      hit->setVolumeID(id);
      // Put it in the collection of hits
      m_intHit->push_back(hit);
      // Keep track of already used id
      m_detectorList[id] = hit;
    }

  // this rotates the hit to local coordinates with respect to the center  
  HepRotation local(*(theTouchable->GetRotation()));
  HepPoint3D center=theTouchable->GetTranslation();

  prePos = local * (prePos-center);
  postPos = local * (postPos-center);
  
  // fill the energy and position    
  hit->addEnergyItem(edep, 
                     McParticleManager::getPointer()->getLastParticle(),
                     (prePos+postPos)/2);
  return true;
}

void IntDetectorManager::EndOfEvent(G4HCofThisEvent*)
{
  // Purpose and Method: this method finalize the manager by saving the
  // integrating hits collection in the TDS
  // Inputs: again G4HCofThisEvent is used only internally by Geant4 and is of
  // no use for us
  // TDS Outputs: The collection of McIntegratingHit is saved in the 
  // /Event/MC/IntegratingHitsCol folder
  
  // Let's sort the hits by volume identifier
  std::sort(m_intHit->begin(),m_intHit->end(), Event::CompareIntHits());

  // store the hits in the TDS
  m_esv->registerObject(EventModel::MC::McIntegratingHitCol , m_intHit);    
#if 0
  // This message is for debug purpouses .. it should be eliminated or converted
  // to a true GAUDI log message
  std::cout << "Actual Event done! " << m_intHit->size() 
            << " integrating hits stored in the TDS" << std::endl;
#endif
}


// File and Version Information:
// $Header$
//
// Description: This is a concrete implementation of the DetectorManager
// abstract class; this one is used to manage sensitive detectors of integrating
// type
//
//
// Author(s):
//      R.Giannitrapani


#include "PosDetectorManager.h"
#include "McParticleManager.h"

#include <iostream>
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "idents/VolumeIdentifier.h"

// Geant4 interface
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"

#include <algorithm>



PosDetectorManager::PosDetectorManager(DetectorConstruction *det,
                                       IDataProviderSvc* esv, IGlastDetSvc* gsv)
  :DetectorManager(det->idMap(), esv, gsv, "PositionDetectorManager")
{
  // See the father class DetectorManager
}

void PosDetectorManager::Initialize(G4HCofThisEvent*)
{
  // Purpose and Method: initialization of the detectors manager
  // Inputs: the G4HCofThisEvent is hinerited from the Geant4 structure and is
  // of no use for our actual implementation of G4Generator (but must be there
  // for Geant4 internal working) 

  // get thr prefix from the GlastDetSvc; if the topvolume is equal to LAT this
  // prefix is empty
  m_prefix = m_gsv->getIDPrefix();

  // clear the list of hited detectors
  m_detectorList.clear();
  // At the start of the event we create a new container
  m_posHit = new Event::McPositionHitVector;    
}

G4bool PosDetectorManager::ProcessHits(G4Step* aStep,
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
  HepPoint3D prePos = aStep->GetPreStepPoint()->GetPosition();
  HepPoint3D postPos = aStep->GetPostStepPoint()->GetPosition();
    
  // determine the ID by studying the history, then call appropriate 
  idents::VolumeIdentifier id = constructId(aStep);

  // If the topvolume is not LAT, prepend a prefix in order to obtain a complete ID
  if (m_prefix.size() != 0)
    id.prepend(m_prefix);

  // Filling of the hits container
  Event::McPositionHit *hit = new Event::McPositionHit;

  // this rotates the hits to local coordinates with respect to the center  
  HepRotation local(*(theTouchable->GetRotation()));
  HepPoint3D center=theTouchable->GetTranslation();

  // this is the global transformation from world to the topVolume; it is the
  // identity if topVolume is equal to LAT
  HepTransform3D global = m_gsv->getTransform3DPrefix();

  McParticleManager* partMan =  McParticleManager::getPointer();
  
  hit->init(edep, id, local*(prePos-center), local*(postPos-center), global*prePos, global*postPos );

  partMan->getLastParticle()->addStatusFlag(Event::McParticle::POSHIT);

  // Track energy at this point
  G4double trkEnergy = aStep->GetTrack()->GetTotalEnergy();
  hit->setParticleEnergy(trkEnergy);

  // Retrieve the particle causing the hit and the ancestor and set the corresponding
  // attributes of the McPositionHit
  if (partMan->getMode() == 1)
    {
      hit->setMcParticle(partMan->getLastParticle());
      hit->setOriginMcParticle(partMan->getOriginParticle());
    }

  // Get the proper time for particle at this hit
  G4double properTime = aStep->GetTrack()->GetProperTime();
  G4double localTime  = aStep->GetTrack()->GetLocalTime();
  hit->setTimeOfFlight(properTime);

  hit->setMcParticleId(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
  hit->setOriginMcParticleId(partMan->getOriginParticle()->particleProperty());

  m_posHit->push_back(hit);
  return true;
}

void PosDetectorManager::EndOfEvent(G4HCofThisEvent*)
{
  // Purpose and Method: this method finalize the manager by saving the
  // integrating hits collection in the TDS
  // Inputs: again G4HCofThisEvent is used only internally by Geant4 and is of
  // no use for us
  // TDS Outputs: The collection of McPositionHit is saved in the 
  // /Event/MC/PositionHitsCol folder
  
  // Let's sort the hits
///////  std::sort(m_posHit->begin(),m_posHit->end(), ComparePosHits());

  // store the hits in the TDS
  m_esv->registerObject( EventModel::MC::McPositionHitCol  , m_posHit);    
#if 0
  // This message is for debug purpouses .. it should be eliminated or converted
  // to a true GAUDI log message
  std::cout << "Actual Event done! " << m_posHit->size() 
            << " position hits stored in the TDS" << std::endl;
#endif
}


#include "McParticleManager.h"
#include "TrackingAction.h"

#include "GlastEvent/MonteCarlo/McParticle.h"

//geant4
#include "G4TrackingManager.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"

//clhep
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"


TrackingAction::TrackingAction()
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  mc::McParticle* parent;
  mc::McParticle* particle;
  McParticleManager* man = McParticleManager::getPointer();
  

#if 0
  G4cout <<  "PRE -- " << aTrack->GetTrackID() << " ";
  G4cout <<  aTrack->GetCurrentStepNumber() << " ";
  G4cout <<  G4BestUnit(aTrack->GetKineticEnergy(), "Energy") << 
    aTrack->GetParentID() << G4endl;
#endif
  
  if (aTrack->GetTrackID() == 1)
    parent = man->getMcParticle(0);
  else
    parent = man->getMcParticle(aTrack->GetParentID());

  particle = new mc::McParticle();
  
  // get the 4-momentum
  
  HepLorentzVector pin(aTrack->GetTotalEnergy(), aTrack->GetMomentum());  
  particle->initialize(parent, aTrack->GetDefinition()->GetPDGEncoding(),
                       mc::McParticle::Swum,pin);
  
  man->addMcParticle(aTrack->GetTrackID(),particle);
  
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  mc::McParticle* particle;
  McParticleManager* man = McParticleManager::getPointer();

#if 0
  G4cout <<  "POS -- " << aTrack->GetTrackID() << " ";
  G4cout <<  aTrack->GetCurrentStepNumber() << " ";
  G4cout <<  G4BestUnit(aTrack->GetKineticEnergy(), "Energy") <<
    G4endl;
#endif

  HepLorentzVector pfin(aTrack->GetTotalEnergy(), aTrack->GetMomentum());  
  particle = man->getMcParticle(aTrack->GetTrackID());
  particle->finalize(pfin, aTrack->GetPosition());
}




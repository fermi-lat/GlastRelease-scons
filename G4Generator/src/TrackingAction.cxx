// File and Version Information:
// $Header$
//
// Description: this method is used to generate new McParticle objects in the
// McParticle hierarchy. It uses a standard mechanism of Geant4 that permits to
// do something whenever a new track is created or an existing track is
// destroyed. To identify the track Geant4 assign an unsigned integer as an
// identifier, starting from 1 for the primary particle
//      
// Author(s):
//      R.Giannitrapani

#include "McParticleManager.h"
#include "TrackingAction.h"

#include "Event/MonteCarlo/McParticle.h"

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
  // Purpose and Method: this method is called every time a new track is created
  // during the Geant4 event simulation
  // Inputs: the G4Track pointer aTrack that gives access to the actual created
  // track object

  Event::McParticle* parent;
  Event::McParticle* particle;
  
  // we get the pointer to the McParticleManager singleton
  McParticleManager* man = McParticleManager::getPointer();

  // if this the primary it has no parent
  if (aTrack->GetTrackID() == 1)
    parent = man->getMcParticle(0);
  else // otherwise it has
    parent = man->getMcParticle(aTrack->GetParentID());

  // lets create a new particle (we don't need to destroy this since it will go
  // in the TDS
  particle = new Event::McParticle();
  
  // get the 4-momentum  
  HepLorentzVector pin(aTrack->GetTotalEnergy(), aTrack->GetMomentum());  
  // we initialize the particle by giving the parent, the PDG encoding, a flag
  // (in that case Swum, and the initial momentum of the particle
  particle->initialize(parent, aTrack->GetDefinition()->GetPDGEncoding(),
                       Event::McParticle::Swum,pin);
  
  // we add this particle to our collection for subsequent saving in the TDS
  man->addMcParticle(aTrack->GetTrackID(),particle);  
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // Purpose and Method: this method is called every time a track is destroied
  // during the Geant4 event simulation
  // Inputs: the G4Track pointer aTrack that gives access to the actual track
  // object
  
  Event::McParticle* particle;

  // we get the pointer to the McParticleManager singleton
  McParticleManager* man = McParticleManager::getPointer();
  // get the 4-momentum   
  HepLorentzVector pfin(aTrack->GetTotalEnergy(), aTrack->GetMomentum());  
  // retrive the particle from our collection with the singleton manager
  particle = man->getMcParticle(aTrack->GetTrackID());
  // we finalize the particle by giving the final momentum and position
  particle->finalize(pfin, aTrack->GetPosition());
}




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
#include "G4Generator/IG4GeometrySvc.h"

#include "idents/VolumeIdentifier.h"

#include "Event/MonteCarlo/McParticle.h"

//geant4
#include "G4TrackingManager.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4TouchableHistory.hh"

//clhep
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include <typeinfo.h>


TrackingAction::TrackingAction(IG4GeometrySvc* gsv):m_geoSvc(gsv)
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Purpose and Method: this method is called every time a new track is created
  // during the Geant4 event simulation
  // Inputs: the G4Track pointer aTrack that gives access to the actual created
  // track object

  bool save = 1;
  Event::McParticle* parent = 0;
  Event::McParticle* particle;
  
  // we get the pointer to the McParticleManager singleton
  McParticleManager* man = McParticleManager::getPointer();
 
  // the name of the process producing the track; by defauld is primary
  std::string process = "primary";

  // if this the primary it has no parent
  if (aTrack->GetTrackID() == 1)
    {
      parent = man->getMcParticle(0);
    }
  else // otherwise it has
    {
	  parent = man->getMcParticle(aTrack->GetParentID());
      process = aTrack->GetCreatorProcess()->GetProcessName();
	}
  
  if ( man->getMode()== 0)
    if ((aTrack->GetTrackID() == 1) || 
        (aTrack->GetParentID() == 1))
      save = 1;
    else save = 0;


  if (save)
    {
      // lets create a new particle (we don't need to destroy this since it will go
      // in the TDS
      particle = new Event::McParticle();
      
      // get the 4-momentum  
      HepLorentzVector pin(aTrack->GetTotalEnergy(), aTrack->GetMomentum());  
      
      
      // we initialize the particle by giving the parent, the PDG encoding, a flag
      // (in that case Swum, and the initial momentum of the particle
      // note that for ions, the stored pdg value: we recreate it from the baryon number 
      G4ParticleDefinition* def = aTrack->GetDefinition();
      int pdgid = def->GetPDGEncoding();
      if ( pdgid==0 ) pdgid = 80000+def->GetBaryonNumber();
      particle->initialize(parent, pdgid,
                           Event::McParticle::Swum,pin,aTrack->GetPosition(),
                           process);
      
      
      idents::VolumeIdentifier ret;
      G4TouchableHistory* theTouchable = (G4TouchableHistory*)aTrack->GetTouchable();
      if (theTouchable)
	  {
	    for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) {
	      const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
//**	      if( physVol->GetMother()==0) break;
	      idents::VolumeIdentifier id = m_geoSvc->getVolumeIdent(physVol);
	      ret.prepend(id);
	    }
	  }
      
      particle->setInitialId(ret);
      
      // we add this particle to our collection for subsequent saving in the TDS
      man->addMcParticle(aTrack->GetTrackID(),particle);  
    }

  // if the particle is an e+ or an e- coming from the conversion of a gamma,
  // than set it as the origin particle, otherwise the primary is the origin
  // particle
  if ((parent == man->getMcParticle(1)) && 
      (parent->particleProperty() == 22) && 
      (process == "conv"))
    man->setOriginParticle(particle);
  else man->setOriginParticle(man->getMcParticle(1));

}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // Purpose and Method: this method is called every time a track is destroied
  // during the Geant4 event simulation
  // Inputs: the G4Track pointer aTrack that gives access to the actual track
  // object


  Event::McParticle* particle = 0;

  // we get the pointer to the McParticleManager singleton
  McParticleManager* man = McParticleManager::getPointer();

  // retrive the particle from our collection with the singleton manager
  particle = man->getMcParticle(aTrack->GetTrackID());

  if (particle)
    {
      // get the 4-momentum   
      HepLorentzVector pfin(aTrack->GetTotalEnergy(), aTrack->GetMomentum()); 

      // Do this for gammas...
      if (particle->particleProperty() == 22 && aTrack->GetTrackStatus() == fStopAndKill)
      {
        // Try getting the final gamma momentum...
        const G4Step*      finalStep          = aTrack->GetStep();
        const G4StepPoint* preFinalStepPoint  = finalStep->GetPreStepPoint();

        pfin = HepLorentzVector(preFinalStepPoint->GetTotalEnergy(), preFinalStepPoint->GetMomentum());
      }

      // we finalize the particle by giving the final momentum and position
      particle->finalize(pfin, aTrack->GetPosition());

      idents::VolumeIdentifier ret;
      G4TouchableHistory* theTouchable = (G4TouchableHistory*)aTrack->GetTouchable();
      if (theTouchable)
	  {
	    for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) {
	      const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
//**	      if( physVol->GetMother()==0) break;
	      idents::VolumeIdentifier id = m_geoSvc->getVolumeIdent(physVol);
	      ret.prepend(id);
	    }
	  }
      
      particle->setFinalId(ret);
      
    }
  
}




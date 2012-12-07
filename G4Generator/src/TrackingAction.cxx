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
#include "McTrajectoryManager.h"
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


TrackingAction::TrackingAction(IG4GeometrySvc* gsv) :
    m_geoSvc(gsv)
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    // Purpose and Method: this method is called every time a new track is created
    // during the Geant4 event simulation
    // Inputs: the G4Track pointer aTrack that gives access to the actual created
    // track object
  
    // we get the pointer to the McParticleManager singleton
    McParticleManager* man = McParticleManager::getPointer();

    if (man->makeMcParticle(aTrack))
    { 
        bool isPrimary = true;

        // the name of the process producing the track; by defauld is primary
        std::string process = "primary";

        // Search for particle, parentage and process info
        Event::McParticle* parent   = man->getMcParticle(aTrack->GetParentID()); // Non-zero unless pruning
        Event::McParticle* particle = man->getMcParticle(aTrack->GetTrackID());  // Is probably zero

        if (aTrack->GetCreatorProcess() != 0) 
        {
            process = aTrack->GetCreatorProcess()->GetProcessName();
            isPrimary = false;
        }

        // lets create a new particle (we don't need to destroy this since it will go
        // in the TDS
        if (particle == 0)
        {
            particle = new Event::McParticle();
      
            // get the 4-momentum  
            CLHEP::HepLorentzVector pin(aTrack->GetTotalEnergy(), aTrack->GetMomentum());  

            // Make sure we have a non-zero parent... this could happen for first particle in list
            if (parent == 0) parent = particle;
      
            // we initialize the particle by giving the parent, the PDG encoding, a flag
            // (in that case Swum, and the initial momentum of the particle
            // note that for ions, the stored pdg value: we recreate it from the baryon number 
            G4ParticleDefinition* def = aTrack->GetDefinition();
            int pdgid = def->GetPDGEncoding();
            if ( pdgid==0 ) pdgid = 80000+def->GetBaryonNumber();
            particle->initialize(parent, pdgid, Event::McParticle::Swum, pin, aTrack->GetPosition(), process);

            // Make sure the primary bit is getting set on the primary particle which is tracked
            if (isPrimary) particle->addStatusFlag(Event::McParticle::PRIMARY);
        }
        // For special runs it can be that the particle already exists
        else
        {
            particle->addStatusFlag(Event::McParticle::Swum);
        }

        idents::VolumeIdentifier ret;
        G4TouchableHistory* theTouchable = (G4TouchableHistory*)aTrack->GetTouchable();
        if (theTouchable)
        {
            for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) 
            {
                const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
                idents::VolumeIdentifier id = m_geoSvc->getVolumeIdent(physVol);
                ret.prepend(id);
            }
        }
      
        particle->setInitialId(ret);
        const G4Material* material = theTouchable->GetVolume(0)->GetLogicalVolume()->GetMaterial();
        particle->setInitialMaterial(material->GetName(), material->GetIndex());
      
        // we add this particle to our collection for subsequent saving in the TDS
        man->addMcParticle(aTrack->GetTrackID(),particle);

        // Now create the McTrajectory to pair with this McParticle
        Event::McTrajectory* trajectory = 
            McTrajectoryManager::getPointer()->getMcTrajectory(aTrack->GetTrackID());
        
        if (!trajectory)
        {
            // A new trajectory to go with this McParticle
            trajectory = new Event::McTrajectory();

            // Add to the collection
            McTrajectoryManager::getPointer()->addMcTrajectory(aTrack->GetTrackID(), trajectory, particle);

            // Make a McTrajectoryPoint corresponding the first point on the track
            float             energy     = aTrack->GetTotalEnergy();
            G4ThreeVector     g4Position = aTrack->GetPosition();
            CLHEP::Hep3Vector point(g4Position);

            Event::McTrajectoryPoint* trajectoryHit = new Event::McTrajectoryPoint(ret, energy, point);

            // Save it
            trajectory->addPoint(trajectoryHit);
        }

        // if the particle is an e+ or an e- coming from the conversion of a gamma,
        // than set it as the origin particle, otherwise the primary is the origin
        // particle
        // Note that this ASSUMES that THE gamma is particle #1 in the McParticle - ID
        // map. I wonder if there is a way to make this work without that assumption?
        if ((parent == man->getMcParticle(1)) && 
            (parent->particleProperty() == 22) && 
            (process == "conv"))
                man->setOriginParticle(particle);
        else man->setOriginParticle(man->getMcParticle(1));
    }
    //else man->setLastParticle(0);
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
    particle = man->getLastParticle();

    if (particle)
    {
        // Should we keep this McParticle?
        if (man->keepMcParticle(aTrack))
        {
            // get the 4-momentum   
            CLHEP::HepLorentzVector pfin(aTrack->GetTotalEnergy(), aTrack->GetMomentum()); 

            // Do this for gammas...
            if (particle->particleProperty() == 22 && aTrack->GetTrackStatus() == fStopAndKill)
            {
                // Try getting the final gamma momentum...
                const G4Step*      finalStep          = aTrack->GetStep();
                const G4StepPoint* preFinalStepPoint  = finalStep->GetPreStepPoint();

                pfin = CLHEP::HepLorentzVector(preFinalStepPoint->GetTotalEnergy(), preFinalStepPoint->GetMomentum());
            }

            // we finalize the particle by giving the final momentum and position
            particle->finalize(pfin, aTrack->GetPosition());

            idents::VolumeIdentifier ret;
            G4TouchableHistory* theTouchable = (G4TouchableHistory*)aTrack->GetTouchable();
            if (theTouchable)
            {
                for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) 
                {
                    const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
                    idents::VolumeIdentifier id = m_geoSvc->getVolumeIdent(physVol);
                    ret.prepend(id);
                }
            }
      
            particle->setFinalId(ret);
            const G4Material* material = theTouchable->GetVolume(0)->GetLogicalVolume()->GetMaterial();
            particle->setFinalMaterial(material->GetName(), material->GetIndex());

            // Move the McParticles to the TDS
            man->saveMcParticle();
            // Also save the trajectories
            McTrajectoryManager::getPointer()->saveMcTrajectory();
        }
        // Otherwise, drop it
        else
        {
            man->dropMcParticle(aTrack);
            McTrajectoryManager::getPointer()->dropMcTrajectory(aTrack);
        }
    }

    return;
}

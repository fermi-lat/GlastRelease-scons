// File and Version Information:
// $Header$
//
// Description: this method is used to check the state of G4 at every step to avoid infinite event loop
// It checks wheter G4 has gone into a G4StateAbort and throw an exception
//      
// Author(s):
//      F.Longo

#include "SteppingAction.h"
#include "G4StateManager.hh"
#include "G4Generator/IG4GeometrySvc.h"
#include "McTrajectoryManager.h"
#include <stdexcept>
#include <string>
#include <sstream>

SteppingAction::SteppingAction(IG4GeometrySvc* geoSvc) : m_geoSvc(geoSvc)
{;}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  if(G4StateManager::GetStateManager()->GetCurrentState()==G4State_Abort) 
    {
      // G4 is going into a G4State Abort
      std::stringstream errorStream;
      errorStream << "G4 in G4StateAbort " ;

      // need to set the Base G4State before aborting the Event Loop

      G4StateManager::GetStateManager()->SetNewState(G4State_Idle);

      throw std::domain_error(errorStream.str());
    }
  // Attempt to add the current point to the McTrajectory (if it exists)
  else
    {
        G4Track*     track = aStep->GetTrack();
        unsigned int id    = track->GetTrackID();

        Event::McTrajectory* trajectory = McTrajectoryManager::getPointer()->getMcTrajectory(id);

        if (trajectory)
        {
            G4ThreeVector g4Position = aStep->GetPostStepPoint()->GetPosition();
            CLHEP::Hep3Vector point(g4Position);

            // Translate the Geant4 volume of this hit into a Glast volume identifier
            idents::VolumeIdentifier volumeID;
            G4TouchableHistory* theTouchable = (G4TouchableHistory*)track->GetTouchable();
            if (theTouchable)
            {
                for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) 
                {
                    const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
                    idents::VolumeIdentifier id = m_geoSvc->getVolumeIdent(physVol);
                    volumeID.prepend(id);
                }
            }

            // Retrieve the total energy at this point
            float energy = track->GetTotalEnergy();

            // Make a McTrajectoryPoint
            Event::McTrajectoryPoint* trajectoryHit = new Event::McTrajectoryPoint(volumeID, energy, point);

            // Save it
            trajectory->addPoint(trajectoryHit);

            // Now deal with creating a relation to either an McPositionHit or McIntegratingHit
            // Retrieve pointer to associated McPositionHit object (if one)
            const Event::McPositionHit* mcPosHit = McTrajectoryManager::getPointer()->getMcPosHit();

            if (mcPosHit)
            {
                Event::McPointToPosHitRel* mcPointHitRel = 
                    new Event::McPointToPosHitRel(trajectoryHit, const_cast<Event::McPositionHit*>(mcPosHit));

                McTrajectoryManager::getPointer()->addMcPosRel(mcPointHitRel);

                // Reset the pointer now that we have used it
                McTrajectoryManager::getPointer()->setMcPosHit(0);
            }

            // Retrieve pointer to associated McPositionHit object (if one)
            const Event::McIntegratingHit* mcIntHit = McTrajectoryManager::getPointer()->getMcIntHit();

            if (mcIntHit)
            {
                Event::McPointToIntHitRel* mcPointHitRel = 
                    new Event::McPointToIntHitRel(trajectoryHit, const_cast<Event::McIntegratingHit*>(mcIntHit));

                McTrajectoryManager::getPointer()->addMcIntRel(mcPointHitRel);

                // Reset the pointer now that we have used it
                McTrajectoryManager::getPointer()->setMcIntHit(0);
            }
        }
    }
}

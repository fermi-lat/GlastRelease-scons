#include "EventAction.hh"
#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"

extern std::ofstream outFile;

TrackingAction::TrackingAction(EventAction* myEA):
eventAction(myEA)
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  if (aTrack->GetTrackID() == 1)
    {
      eventAction->SetStartEn(aTrack->GetKineticEnergy());
      eventAction->SetInitialDir(aTrack->GetMomentumDirection());
    }
  
#if 0
  if ((aTrack->GetDefinition()->GetParticleName() == "e-") && 
      (aTrack->GetParentID() == 1) &&
      (aTrack->GetCreatorProcess()->GetProcessName() == "conv"))
    {
      outFile << "1 " << (aTrack->GetKineticEnergy())/MeV << G4endl;
    }

  if ((aTrack->GetDefinition()->GetParticleName() == "e+") && 
      (aTrack->GetParentID() == 1) &&
      (aTrack->GetCreatorProcess()->GetProcessName() == "conv"))
    outFile << "2 " << (aTrack->GetKineticEnergy())/MeV << G4endl;
#endif
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  if (aTrack->GetTrackID() == 1)
    {
      eventAction->SetEndEn(aTrack->GetKineticEnergy());
      eventAction->SetFinalDir(aTrack->GetMomentumDirection());
    }

}




#include "TrackingAction.h"
#include "G4TrackingManager.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"

TrackingAction::TrackingAction()
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
#if 0
  G4cout <<  "PRE -- " << aTrack->GetTrackID() << " ";
    G4cout <<  aTrack->GetCurrentStepNumber() << " ";
    G4cout <<  G4BestUnit(aTrack->GetKineticEnergy(), "Energy") << 
       G4endl;
#endif
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
#if 0
  G4cout <<  "POS -- " << aTrack->GetTrackID() << " ";
    G4cout <<  aTrack->GetCurrentStepNumber() << " ";
    G4cout <<  G4BestUnit(aTrack->GetKineticEnergy(), "Energy") <<
     G4endl;
#endif
}




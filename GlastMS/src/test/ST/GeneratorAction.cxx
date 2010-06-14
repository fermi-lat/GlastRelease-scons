// --------------------------------------------------------------
//      GEANT 4 -  2000
//
// Comments
//
// - implementation file of event manager
//
//---------------------------------------------------------------

#include "GeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"

GeneratorAction::GeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/gun/particle e+");
  UI->ApplyCommand("/gun/energy 1.0 GeV");
  UI->ApplyCommand("/gun/position -0.1 0.0 0.0 m");
  UI->ApplyCommand("/gun/direction 1.0 0.0 0.0");
}

GeneratorAction::~GeneratorAction()
{
  delete particleGun;
}

void GeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int i = anEvent->GetEventID();
  particleGun->GeneratePrimaryVertex(anEvent);
}








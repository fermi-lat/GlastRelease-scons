// --------------------------------------------------------------
//      GEANT 4 -  2000
//
// by Riccardo Giannitrapani (august 2000)
// Credits to Francesco Longo and to the Geant4 user support;
// parts of this simulation is taken by standard examples in 
// the Geant4 distribution and from tb2000 simulation by F.Longo  
// 
// Credits to Glast ground software team for suggestions and
// help
// 
// Comments
//
// - code to test GEANT4 physical processes in a very simple
//   simulation; just a particle scattering a box of a certain
//   material
// - Hits are stored in a very hugly way, by defining a global
//   file stream and putting data in it. The data saved are 
//   energy deposition in the sensible material (this come from the 
//   LayerHit class) and energy loss of the particle (this come from
//   TrackingAction, but are stored in the LayerHit class, that's 
//   not a good thing, so I have to change it). Units are MeV. 
//   This file dumping has to change in the future or at least 
//   has to be done better.
// - Since this simulation has been built in different times and
//   places, I've not used a very good naming convention. For
//   example the hit class is called LayerHit, since in a first
//   moment I was working on a layred geometry. Again this has
//   to change in the future.
// - The macro file prerun.mac is loaded by default; you can
//   put your values for particle type, energy, material type, 
//   thickness of the material, number of runs and so on in this file
// ----------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "VisManager.hh"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "GeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"

// This global file is used to store relevant data for
// analysis with a separate program
std::ofstream outFile("test.dat");

int main()
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);

  G4UIsession* session=0;
  session = new G4UIterminal();

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new VisManager;
  visManager->Initialize();
#endif

  // set mandatory user action classes
  runManager->SetUserAction(new GeneratorAction);
  runManager->SetUserAction(new RunAction);

  // set optional user action classes
  EventAction* eventAction = new EventAction;
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new TrackingAction(eventAction));

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      // prerun.mac is loaded by default
      UI->ApplyCommand("/control/execute src/test/ST/prerun.mac");    
      session->SessionStart();
      delete session;
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  return 0;
}



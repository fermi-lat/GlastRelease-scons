#include "EventAction.hh"
#include "LayerHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

// This file is a global variabe in which I store energy deposition per hit
// It is a sort of hack, I think I have to use some other way to export 
// hit information ....
extern std::ofstream outFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction()
  :printModulo(10),layerCollID(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void EventAction::BeginOfEventAction(const G4Event* evt)
{
  
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0)
   { 
     G4cout << "\n---> Event: " << evtNb << G4endl;
   }
    
 if (layerCollID==-1)
 {
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  layerCollID = SDman->GetCollectionID("LayerCollection");
 } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();


  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  LayerHitsCollection* CHC = NULL;
  if (HCE)
      CHC = (LayerHitsCollection*)(HCE->GetHC(layerCollID));
       
    if (CHC)
      {
	int n_hit = CHC->entries();
	G4double totESil=0, totLSil=0;
	for (int i=0;i<n_hit;i++)
	  {
	    totESil += (*CHC)[i]->GetEdepSil(); 
	    totLSil += (*CHC)[i]->GetTrakSil();
	  }
	
	// Here I put the total energy deposition in the file for 
	// later analysis
        outFile << totESil/MeV << " " << (startEn - endEn)/MeV 
                << " " << finalDir.angle(initialDir)/deg << G4endl;
      }
    
}








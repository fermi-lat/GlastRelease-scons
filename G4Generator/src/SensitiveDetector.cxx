#include <map>

#include "DetectorConstruction.h"
#include "SensitiveDetector.h"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SensitiveDetector::SensitiveDetector(G4String name,
                                   DetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SensitiveDetector::~SensitiveDetector()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool SensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  typedef map<G4VPhysicalVolume*, string> M;
  M::const_iterator i; 

  // Energy Deposition & Step Lenght

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double stepl = aStep->GetStepLength();

  if ((edep==0.)) return false;      

  // Physical Volume
  
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4LogicalVolume* logVol = physVol->GetLogicalVolume();
  G4String material = logVol->GetMaterial()->GetName();
  G4String nameVolume = physVol->GetName();
  G4String name = "";

  do
    { 
      i = Detector->g4ID.find(physVol);
      if(i != Detector->g4ID.end()){
	name = i->second + name;
      }
      theTouchable->MoveUpHistory();
      physVol = theTouchable->GetVolume(); 
    } while(physVol->GetName() != "motherVolume");

  
  // Initial Position

  G4ThreeVector InitPos = aStep->GetPreStepPoint()->GetPosition();
  
  // Final Position 
  
  G4ThreeVector FinPos = aStep->GetPostStepPoint()->GetPosition();

  std::cout << "Hit -> " << " " << 
    aStep->GetTrack()->GetDefinition()->GetParticleName() << 
    " " << 
    edep << 
    " " << 
    stepl << 
    " " <<
    name << " " <<
    nameVolume << " " << material << " " << 
    std::endl;

  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{

}






//$Header$

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
    m_energySummary.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool SensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
    typedef std::map<G4VPhysicalVolume*, std::string> M;
  M::const_iterator i; 

  // Energy Deposition & Step Lenght

  G4double edep = aStep->GetTotalEnergyDeposit()/MeV;
  G4double stepl = aStep->GetStepLength()/mm;

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
#if 1
  std::cout << "Hit -> " 
	    << std::setw(8) << aStep->GetTrack()->GetDefinition()->GetParticleName()  
	    << std::setw(8) << std::setprecision(3) << edep  
	    << std::setw(8) << std::setprecision(3) <<stepl 
	    << std::setw(18) << name << " " 
	    << std::setw(12) << nameVolume << " " << material
	    << std::endl;
#endif
  m_energySummary[logVol->GetName()] += edep;
  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{

    std::cout << "Energy deposit summary: " << std::endl;
    for( std::map<std::string, double>::const_iterator it = m_energySummary.begin(); it!= m_energySummary.end(); ++it){
        std::cout 
            << std::setw(20) << (*it).first 
            << std::setw(10) << std::setprecision(3) << (*it).second << std::endl;
    }

}






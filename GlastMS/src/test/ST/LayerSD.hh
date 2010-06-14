#ifndef LayerSD_h
#define LayerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "LayerHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class LayerSD : public G4VSensitiveDetector
{
  public:
  
      LayerSD(G4String, DetectorConstruction* );
     ~LayerSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      LayerHitsCollection*  LayerCollection;      
      DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

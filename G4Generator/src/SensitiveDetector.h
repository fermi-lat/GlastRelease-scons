#ifndef SensitiveDetector_h
#define SensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"


class DetectorConstruction;
class G4HCofThisEvent;
class G4Step;

class SensitiveDetector : public G4VSensitiveDetector
{
 public:
  
  SensitiveDetector(G4String, DetectorConstruction* );
  ~SensitiveDetector();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear(){};
  void DrawAll(){};
  void PrintAll(){};
  
 private:
  DetectorConstruction* Detector;
};

#endif

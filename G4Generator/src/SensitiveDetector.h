//$Header$

#ifndef SensitiveDetector_h
#define SensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include <map>
#include <string>

class DetectorConstruction;
class G4HCofThisEvent;
class G4Step;

class SensitiveDetector : public G4VSensitiveDetector
{
 public:
  
  SensitiveDetector(G4String, DetectorConstruction* );
  ~SensitiveDetector();
  
  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);
  //void clear(){};
  
 private:
  DetectorConstruction* Detector;

  // for debugging: summary of energy per logical volume
  std::map<std::string, double> m_energySummary;
};

#endif

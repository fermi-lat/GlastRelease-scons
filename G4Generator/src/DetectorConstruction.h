#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4VPhysicalVolume;
class SensitiveDetector;

#include <string>
#include <map>
#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  DetectorConstruction();
  ~DetectorConstruction();
  
 public:
  G4VPhysicalVolume* Construct();

  std::map <G4VPhysicalVolume*, std::string> g4ID;    

 private:
  void DefineMaterial();
  G4VPhysicalVolume* GeometryConstruct();
  
  SensitiveDetector* sDetector;  
};

#endif


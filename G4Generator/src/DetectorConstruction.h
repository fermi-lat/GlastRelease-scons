#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

 private:
     void DefineMaterial();
     G4VPhysicalVolume* GeometryConstruct();
};

#endif


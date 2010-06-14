// --------------------------------------------------------------
//      GEANT 4 - GLAST  2000
//---------------------------------------------------------------


#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1


#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;
class LayerSD;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();
  
public:
  G4VPhysicalVolume* Construct();
  void DefineMaterial();
  G4VPhysicalVolume* GeometryConstruct();
  
  void SetMaterial(G4String);
  void SetSizeY(G4double);
  void SetSizeZ(G4double);
  void SetThickness(G4double);
  void UpdateGeometry();
  
private:
  G4Material* Material;
  G4Material* HallMaterial;
  G4double Thickness;
  G4double SizeY;
  G4double SizeZ;
  G4LogicalVolume* Tile_log;


  DetectorMessenger* Messenger; // pointer to the messenger
  LayerSD* layerSD;  //pointer to the sensitive detector of the single layer
  
};

#endif


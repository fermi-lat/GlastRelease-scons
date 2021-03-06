// --------------------------------------------------------------
//      GEANT 4 - GLAST
// --------------------------------------------------------------

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "LayerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "globals.hh"


// --------- CONSTRUCTOR -----------------------------------
DetectorConstruction::DetectorConstruction()
:layerSD(NULL),Material(NULL),HallMaterial(NULL),Tile_log(NULL)
{
  Thickness = 2.35*cm;
  SizeY = 100*cm;
  SizeZ = 100*cm;

  Messenger = new DetectorMessenger(this);
}

// --------- DESTRUCTOR -------------------------------------

DetectorConstruction::~DetectorConstruction()
{;}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterial();
  return GeometryConstruct();
}


void DetectorConstruction::DefineMaterial()
{
  //------- MATERIALS --------------------------------------

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4String name, symbol;
  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  G4double temperature, pressure;


  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Si", z=14., a, density);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 12.011*g/mole;
  G4Element* Ca  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 1.008*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen"  ,symbol="H" , z= 1., a);

  a = 26.982*g/mole;
  G4Element* Alu  = new G4Element(name="Aluminum"  ,symbol="Al" , z= 13., a);

  density = 1.64*g/cm3;
  G4Material* rubber60A1 = 
    new G4Material(name="rubber60A1", density, ncomponents=3);
  rubber60A1->AddElement(Ca, natoms=1);
  rubber60A1->AddElement(H, natoms=1);
  rubber60A1->AddElement(Alu, natoms=1);


  density = 1.0320*g/cm3;
  G4Material* Polystyrene = 
    new G4Material(name="Polystyrene", density, ncomponents=2);
  Polystyrene->AddElement(Ca, natoms=1);
  Polystyrene->AddElement(H, natoms=1);

  density = 5.0E-02*g/cm3;
  G4Material* FOAM05 = 
    new G4Material(name="FOAM05", density, ncomponents=2);
  FOAM05->AddElement(Ca, natoms=1);
  FOAM05->AddElement(H, natoms=1);


  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Al", z=13., a, density);

  density = 5.4040E-01*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al20 = new G4Material(name="Al20", z=13., a, density);

  density = 1.3520*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al50 = new G4Material(name="Al50", z=13., a, density);

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  density = 2.265*g/cm3;
  a = 12.01*g/mole;
  G4Material* C  = new G4Material(name="C"  , z= 6., a, density);

  density = 19.3*g/cm3;
  a = 183.84*g/mole;
  G4Material* W = new G4Material(name="W", z=74., a, density);


  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;

  G4Material* Galactic = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

  G4Material* Vacuum = new G4Material(name="Vacuum", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);
  // TODO

  a = 132.90545*g/mole;
  G4Element* Cs  = new G4Element(name="Cesium",symbol="Cs" , z= 55., a);

  a = 126.90447*g/mole;
  G4Element* I  = new G4Element(name="Iodine",symbol="I" , z= 53., a);

  density = 4.53*g/cm3;
  G4Material* CsI = 
    new G4Material(name="CsI", density, ncomponents=2);
  CsI->AddElement(Cs, natoms=1);
  CsI->AddElement(I, natoms=1);

}

G4VPhysicalVolume* DetectorConstruction::GeometryConstruct()
{
  //---------- VOLUMES --------------------------------------------
  //------------------------------ beam line along x axis

  //--------- EXPERIMENTAL HALL ---------------------------

  G4Material* ptMaterial = G4Material::GetMaterial("Galactic");


  G4double expHall_x = 1.0*m;
  G4double expHall_y = 1.0*m;
  G4double expHall_z = 1.0*m;

  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* experimentalHall_log
    = new
G4LogicalVolume(experimentalHall_box,ptMaterial,"expHall_log",0,0,0);

  G4VPhysicalVolume* experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall",
                       experimentalHall_log,0,false,0);

  // ---------- MATERIAL TILE -------------------------------

  G4Box* Tile_box
    = new G4Box("Tile_box",Thickness/2,SizeY/2,SizeZ/2);

  Tile_log= new G4LogicalVolume(Tile_box,Material,"Tile_log",0,0,0);

  G4double TilePos_x = 0.0*m;
  G4double TilePos_y = 0.0*m;
  G4double TilePos_z = 0.0*m;

  G4VPhysicalVolume* Tile_phys
    = new G4PVPlacement(0,
			G4ThreeVector(TilePos_x, TilePos_y, TilePos_z),
			Tile_log,"Tile",experimentalHall_log,false,0);
  


  //------------- SENSITIVE VOLUMES ---------------------------------------


  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!layerSD)
  {
    layerSD = new LayerSD("LayerSD",this);
    SDman->AddNewDetector( layerSD );
  }

  if (Tile_log)
      Tile_log->SetSensitiveDetector(layerSD);
   

  return experimentalHall_phys;
}


void DetectorConstruction::SetMaterial(G4String mat)
{
  G4Material* ptMaterial = G4Material::GetMaterial(mat);
  if(ptMaterial)
    {
      Material = ptMaterial;
      Tile_log->SetMaterial(ptMaterial); 
    }
}

void DetectorConstruction::SetThickness(G4double value)
{
  Thickness = value;
}

void DetectorConstruction::SetSizeY(G4double value)
{
  SizeY = value;
}

void DetectorConstruction::SetSizeZ(G4double value)
{
  SizeZ = value;
}

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(GeometryConstruct());
}









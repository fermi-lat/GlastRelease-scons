// $Header$
#include "DetectorConstruction.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SDManager.hh"

#include "GlastDetectorManager.h"

#include "G4Geometry.h"
#include "G4Media.h"

#include <iomanip>
#include <cassert>

DetectorConstruction::DetectorConstruction(IGlastDetSvc* gsv,
                                           IDataProviderSvc* esv):m_gsv(gsv)
{
    
  // now create the GlastDetector manager, and pass in the id map
  m_glastdet = new GlastDetectorManager(this, esv);


}

DetectorConstruction::~DetectorConstruction()

{ 
  delete m_glastdet;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  G4Geometry* geom = new G4Geometry(m_glastdet, &m_idMap);
  G4Media* media = new G4Media();
  m_gsv->accept(*media);
  m_gsv->accept(*geom);
  //  std::cout << (*G4Material::GetMaterialTable()) <<std::endl;
  std::cout << "Geometry done with " << geom->getPhysicalNumber() << 
      " physical volumes" << std::endl;
  return geom->getWorld();
}

















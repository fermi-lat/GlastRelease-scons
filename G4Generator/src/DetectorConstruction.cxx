#include "DetectorConstruction.h"
#include "SensitiveDetector.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SDManager.hh"

#include "detModel/Management/Manager.h"
#include "detModel/Management/XercesBuilder.h"
#include "detModel/Sections/Shape.h"
#include "detModel/Gdd.h"

#include "G4SectionsVisitor.h"
#include "G4MaterialsVisitor.h"
#include "GlastDetectorManager.h"

#include <iomanip>

DetectorConstruction::DetectorConstruction(std::string topvol, std::string visitorMode)
: m_topvol(topvol)
{
  std::string filename= std::string(::getenv("XMLUTILROOT"))+"/xml/flight.xml" ;

  detModel::Manager* gddManager = detModel::Manager::getPointer();
  
  gddManager->setBuilder(new detModel::XercesBuilder);  

  gddManager->setNameFile(filename);

  gddManager->setMode(visitorMode);

  gddManager->build(detModel::Manager::all);
  // now create the GlastDetector manager, to pass in the id map
  m_glastdet = new GlastDetectorManager(this);
  
}

DetectorConstruction::~DetectorConstruction()
{ 
  detModel::Manager* gddManager = detModel::Manager::getPointer();
  gddManager->cleanGdd(); 
  delete m_glastdet;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // This method build the material -- MUST be implemented with a MaterialVisitor
  detModel::Manager::getPointer()->startVisitor(&G4MaterialsVisitor());
  
  // This method build the geometry with a visitor
  return GeometryConstruct();
}


G4VPhysicalVolume* DetectorConstruction::GeometryConstruct()
{
 

  detModel::Manager* gddManager = detModel::Manager::getPointer();
  detModel::Gdd* gdd = gddManager->getGdd();
  
  // create the visitor as with the topvolume, with a pointer to the map to fill
  G4SectionsVisitor visitor (m_topvol, &m_idMap);
  

  gddManager->startVisitor(&visitor);
  
  for( std::vector<G4LogicalVolume*>::iterator it = visitor.g4Logicals.begin(); 
  it != visitor.g4Logicals.end(); ++it){
      G4LogicalVolume * lv = *it;
      detModel::Shape* sh = gdd->getShapeByName(lv->GetName());
      if( sh  && sh->getSensitive() ){
          std::cout << "Found sensitive: " << sh->getName() << std::endl;
          m_glastdet->process(lv);
      }
  }
  //----------------------------------------


   G4VPhysicalVolume* res = visitor.worldphys;
  visitor.summary(std::cout);
  
  return res;
}

















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

#include "detModel/Management/Manager.h"
#include "detModel/Management/XercesBuilder.h"
#include "detModel/Sections/Shape.h"
#include "detModel/Sections/Volume.h"
#include "detModel/Sections/Box.h"
#include "detModel/Gdd.h"
#include "detModel/Management/IDmapBUilder.h"
#include "detModel/Utilities/PositionedVolume.h"

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
  
  detModel::IDmapBuilder idmap(topvol);
  
  gddManager->startVisitor(&idmap);
  idmap.summary(std::cout);

  DisplayManager* dm = DisplayManager::instance();
  if (dm !=0){
      for( detModel::IDmapBuilder::PVmap::const_iterator id = idmap.begin(); id!=idmap.end(); ++id){
          const detModel::PositionedVolume * pv = (*id).second;
          const detModel::Volume* vol = pv->getVolume();
          const detModel::Box* b = dynamic_cast<const detModel::Box*>(pv->getVolume());
          if (b !=0) {
              dm->addDetectorBox(HepTransform3D( pv->getRotation(),pv->getTranslation()), b->getX(), b->getY(), b->getZ());
          }
      }
  }

}

DetectorConstruction::~DetectorConstruction()

{ 
  detModel::Manager* gddManager = detModel::Manager::getPointer();
  gddManager->cleanGdd(); 
  delete m_glastdet;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  detModel::Manager* gddManager = detModel::Manager::getPointer();

  //  build the material  with a MaterialVisitor
  gddManager->startVisitor(&G4MaterialsVisitor());
  
  //  build the geometry with a visitor
  // create the visitor as with the topvolume, with a pointer to the map to fill
  G4SectionsVisitor visitor (m_topvol, &m_idMap);
  
  gddManager->startVisitor(&visitor);

  detModel::Gdd* gdd = gddManager->getGdd();

  for( G4SectionsVisitor::Logicals::const_iterator it = visitor.begin(); it !=visitor.end(); ++it){
      G4LogicalVolume * lv = *it;
      detModel::Shape* sh = gdd->getShapeByName(lv->GetName());
      if( sh  && sh->getSensitive() ){
          m_glastdet->process(lv);
      }
  }
  //----------------------------------------


  //G4VPhysicalVolume* res = visitor.worldphys;
   //visitor.summary(std::cout);
  
  return visitor.getWorld();
}

















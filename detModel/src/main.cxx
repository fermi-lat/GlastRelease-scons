/// Test program for the generic model

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "detModel/Management/Manager.h"
#include "detModel/Management/VrmlSectionsVisitor.h"
#include "detModel/Management/CountMaterial.h"
#include "detModel/Management/IDmapBuilder.h"
#include "detModel/Management/PrinterSectionsVisitor.h"
#include "detModel/Management/HtmlConstantsVisitor.h"
#include "detModel/Management/PrinterMaterialsVisitor.h"
#include "detModel/Management/XercesBuilder.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Materials/Material.h"
#include "detModel/Utilities/Color.h"
#include "detModel/Utilities/PositionedVolume.h"
#include "detModel/Sections/Volume.h"
#include "detModel/Sections/Shape.h"
#include "detModel/Sections/Box.h"
#include "detModel/Gdd.h"


/* This basic test needs two argument; the xml file to use and the volume 
   name to use as the mother volume (for example oneCAL). If only the filename
   is specified, the topVolume specified in the section is used as the mother
   volume.

   Ex.
   ./test.exe ../../../xmlUtil/v2r1/xml/flight.xml oneTKR
   
   The test program produce a sections.wrl file with the VRML representation 
   of the geometry (be careful, if you have a big geometry the file can be
   huge) and a file constants.html with the constants tables in html format. 
*/
int main(int argc, char* argv[]) {

  // We need the XML flight as input to the test executable
  if (argc == 1) return 0;

  // We retrive the manager pointer (it is a singleton, so it is not possible
  // to create it in the usual way)
  detModel::Manager* manager = detModel::Manager::getPointer();

  // Set the builder and the file name
  manager->setBuilder(new detModel::XercesBuilder);
  manager->setNameFile(argv[1]);
  
  // We set the mode for the choice elements in the XML file
  //manager->setMode("digi");

  // We build the hierarchy; in that case we build all, i.e. both the constants
  // the sections and the materials
  manager->build(detModel::Manager::all);

  // We start the VRMLSectionsVisitor to build the vrml file
  // If we don't specify a string in the constructor, it will build the
  // vrml file for all the volumes placed in the topVolume, otherwise it
  // will build the vrml file for the specified volume. The output is
  // placed in sections.vrml

  /*
  detModel::VrmlSectionsVisitor* visitor;
  if (argc == 2)
    visitor = new detModel::VrmlSectionsVisitor("");  
  else
    visitor = new detModel::VrmlSectionsVisitor(argv[2]);  

  visitor->setMode("digi");
  */


  // We retrive the hierarchy entry point, i.e. the GDD object. It
  // contains all the relevant information
  detModel::Gdd* g = manager->getGdd();
  
  // An example; we retrive some info from the xml file
  std::cout << "XML file contains " << g->getVolumesNumber() << " volumes." << std::endl;
  std::cout << "XML file contains " << g->getMaterialsNumber() << " materials." << std::endl;
  std::cout << "XML file contains " << g->getConstantsNumber() << " constants." << std::endl;

  double maxLogVal;
  bool getOk;
  getOk = manager->getNumericConstByName("maxLog", &maxLogVal);
  if (getOk)
    std::cout << "(Float) value for maxLog is " << maxLogVal << std::endl;
  else std::cout << "getNumericConstByName failed" << std::endl;

  int maxLogInt;
  getOk = manager->getNumericConstByName("maxLog", &maxLogInt);
  if (getOk)
    std::cout << "(int) value for maxLog is " << maxLogInt << std::endl;
  else std::cout << "getNumericConstByName (int) failed" << std::endl;

  // -------
  double nCsISeg;
  getOk = manager->getNumericConstByName("nCsISeg", &nCsISeg);
  if (getOk)
    std::cout << "(Float) value for nCsISeg is " << nCsISeg << std::endl;
  else std::cout << "getNumericConstByName failed" << std::endl;

  int nCsISegInt;
  getOk = manager->getNumericConstByName("nCsISeg", &nCsISegInt);
  if (getOk)
    std::cout << "(int) value for nCsISeg is " << nCsISegInt << std::endl;
  else std::cout << "getNumericConstByName (int) failed" << std::endl;

  double xNum;
  getOk = manager->getNumericConstByName("xNum", &xNum);
  if (getOk)
    std::cout << "(Float) value for xNum is " << xNum << std::endl;
  else std::cout << "getNumericConstByName failed" << std::endl;

  int xNumInt;
  getOk = manager->getNumericConstByName("xNum", &xNumInt);
  if (getOk)
    std::cout << "(int) value for xNum is " << xNumInt << std::endl;
  else std::cout << "getNumericConstByName (int) failed" << std::endl;





  // -------

  double diodeX;
  getOk = manager->getNumericConstByName("diodeX", &diodeX);
  if (getOk)
    std::cout << "(Float) value for diodeX is " << diodeX << std::endl;
  else std::cout << "getNumericConstByName failed" << std::endl;

  int diodeXInt;
  getOk = manager->getNumericConstByName("diodeX", &diodeXInt);
  if (getOk)
    std::cout << "(int) value for diodeX is " << diodeXInt << std::endl;
  else std::cout << "getNumericConstByName (int) failed" << std::endl;

  
  // Retrive the materials, generate the colors and set some transparency values
  detModel::MatCollection* mats = g->getMaterials();  
  mats->generateColor();
  mats->setMaterialTransparency("Vacuum",0.7);
  
  // We start the HTMLConstantsVisitor to build the html file with the
  // constants tables. Colors and layout are stolen from Joanne ones.
  manager->startVisitor(new detModel::HtmlConstantsVisitor());

  // We set a mode for choices
  //  manager->setMode("propagate");
  // We start the vrml visitor
  /*
  manager->startVisitor(visitor);
  */
  manager->startVisitor(new detModel::CountMaterial(argv[2]));

  detModel::IDmapBuilder idMap(argv[2]);
  //idMap.setMode("digi");
  manager->startVisitor(&idMap);
  idMap.summary(std::cout);

  // Try out new IDmapBuilder methods

  std::cout << "Top volume is " << std::string(argv[2]) << std::endl;
  idents::VolumeIdentifier worldId;
  idents::VolumeIdentifier relId;

  relId.append(1);     // fTowerObjects field
  relId.append(1);     // tray number
  relId.append(0);     // fmeas = x
  relId.append(1);     // Si top

  idents::VolumeIdentifier prefix = idMap.getIDPrefix();

  std::cout << "ID prefix is " << prefix.name() << std::endl;

  worldId = prefix;
  worldId.append(relId);

  std::string  s;
  std::vector<double> params;
  bool ret = idMap.getShapeByID(worldId, &s, &params);
  if (ret) {
    std::cout << "String from getShapeByID(worldId..) is " << s << std::endl;
  }
  else std::cout << "getShapeByID call failed" << std::endl;

  ret = idMap.getShapeByID(relId, &s, &params);
  if (ret) {
    std::cout << "String from getShapeByID(relId..) is " << s << std::endl;
  }
  else std::cout << "getShapeByID call failed" << std::endl;



  ret = idMap.getShapeByTopID(worldId, &s, &params);

  if (ret) {
    std::cout << "String from getShapeByTopID(worldId..) is " 
              << s << std::endl;
  }
  else std::cout << "getShapeByTopID call failed" << std::endl;

  ret = idMap.getShapeByTopID(relId, &s, &params);

  if (ret) {
    std::cout << "String from getShapeByTopID(relId..) is " 
              << s << std::endl;
  }
  else std::cout << "getShapeByTopID call failed" << std::endl;


  HepTransform3D   trans;
  ret = idMap.getTransform3DByID(worldId, &trans);

  ret = idMap.getTopTransform3DByID(relId, &trans);
  /* Diagnostic */
  //  for( detModel::IDmapBuilder::PVmap::const_iterator id = idMap.begin(); 
  //       id!=idMap.end(); ++id){
  //    const detModel::PositionedVolume * pv = (*id).second;
  //    const detModel::Volume* vol = pv->getVolume();
  //    const detModel::Box* b = 
  //      dynamic_cast<const detModel::Box*>(pv->getVolume());
  //    if (b !=0) {
  //      std::cout << "Got box with id ";
  //       std::cout << ((*id).first).name() << " named "
  //                << vol->getName() << std::endl;
  //    }
  //  }

  /*
  delete visitor;
  */
  delete manager;
  return(0);
}






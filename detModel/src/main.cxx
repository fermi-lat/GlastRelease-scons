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
#include "idents/VolumeIdentifier.h"
#include "idents/AcdId.h"

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

// 
void orderSegments(std::vector<idents::VolumeIdentifier>& segs, 
                   short face, short ribNum, bool xOrient, bool increasing,
                   detModel::IDmapBuilder& idMap);

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

  detModel::VrmlSectionsVisitor* visitor;

  std::string volName("");
  if (argc > 2) volName = std::string(argv[2]);

  visitor = new detModel::VrmlSectionsVisitor(volName);
  visitor->setMode("digi");

  // We retrieve the hierarchy entry point, i.e. the GDD object. It
  // contains all the relevant information
  detModel::Gdd* g = manager->getGdd();
  
  // An example; we retrieve some info from the xml file
  std::cout << "XML file contains " << g->getVolumesNumber() << " volumes." << std::endl;
  std::cout << "XML file contains " << g->getMaterialsNumber() << " materials." << std::endl;
  std::cout << "XML file contains " << g->getConstantsNumber() << " constants." << std::endl;

  double maxLogVal;
  manager->getNumericConstByName("maxLog", &maxLogVal);
  std::cout << "value for maxLog is " << maxLogVal << std::endl;

  // Retrieve the materials, generate the colors and set some 
  // transparency values

  detModel::MatCollection* mats = g->getMaterials();  
  mats->generateColor();
  mats->setMaterialTransparency("Vacuum",0.7);

  // We start the HTMLConstantsVisitor to build the html file with the
  // constants tables. Colors and layout are stolen from Joanne ones.

  manager->startVisitor(new detModel::HtmlConstantsVisitor());
  
  // We set a mode for choices
  //  manager->setMode("propagate");
  // We start the vrml visitor
              /*   Comment out for now
               manager->startVisitor(visitor);
                 */
  manager->startVisitor(new detModel::CountMaterial(volName));

  detModel::IDmapBuilder idMap(volName);
  manager->startVisitor(&idMap);

  idMap.summary(std::cout);
  
  std::vector<idents::VolumeIdentifier> segments;
  orderSegments(segments, 1, 1, true, true, idMap);
  segments.clear();

  orderSegments(segments, 1, 3, false, true, idMap);
  segments.clear();

  orderSegments(segments, 1, 1, true, false, idMap);
  segments.clear();

  orderSegments(segments, 2, 1, true, true, idMap);
  segments.clear();

  orderSegments(segments, 2, 1, false, true, idMap);
  segments.clear();

  orderSegments(segments, 0, 1, true, true, idMap);
  segments.clear();

  orderSegments(segments, 0, 1, false, true, idMap);
  segments.clear();


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

  delete visitor;
  delete manager;
  return(0);
}

// Something like this will be added to GlastSvc
void orderSegments(std::vector<idents::VolumeIdentifier>& segs,
                   short face, short ribNum, bool xOrient, bool increasing,
                   detModel::IDmapBuilder& idMap) {

  using detModel::IDmapBuilder;

  // Make a volume identifier for a ribbon
  idents::VolumeIdentifier sample;
  IDmapBuilder::axisDir dir = IDmapBuilder::zDir; // except for face=0

  sample.append(1);           // ACD
  sample.append(face);
  sample.append(41);   // ribbon

  if (xOrient) {
    sample.append(1);  // measures Y
    if (face==0) dir = IDmapBuilder::xDir;
  }
  else {
    sample.append(0);  // measures X
    if (face==0) dir = IDmapBuilder::yDir;
  }
 
  sample.append(ribNum);
  sample.append(0);    // segment # doesn't matter; will be masked off

  std::vector<bool> wild;
  wild.reserve(6);
  wild.push_back(false);
  wild.push_back(false);
  wild.push_back(false);
  wild.push_back(false);
  wild.push_back(false);
  wild.push_back(true);

  idMap.orderSensitive(segs, sample, wild, dir, increasing);

  std::cout << "For inputs face = " << face << " ribNum = " << ribNum
            << " xOrient = " << xOrient << " and increasing = "
            << increasing << std::endl
            << "output vol ids are: " << std::endl;
  for (unsigned ix = 0; ix < segs.size(); ix++) {
    std::cout <<  segs[ix].name() << std::endl;
  }
  std::cout << std::endl;

}





/// Test program for the generic model

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDVRMLSectionsVisitor.h"
#include "detModel/Management/GDDHTMLConstantsVisitor.h"
#include "detModel/Management/GDDXercesBuilder.h"
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDshape.h"
#include "detModel/Sections/GDDbox.h"
#include "detModel/GDD.h"


int main(int argc, char* argv[]) {

  // We need the XML flight as input to the test executable
  if (argc == 1) return 0;

  // We retrive the manager pointer (it is a singleton, so it is not possible
  // to create it in the usual way)
  GDDmanager* manager = GDDmanager::getPointer();
 
  // We set the builder; the XercesBuilder needs the name of the XML file
  // in the constructor
  manager->setBuilder(new GDDXercesBuilder(argv[1]));

  // We set the mode for the choice elements in the XML file
  manager->setMode("propagate");

  // We build the hierarchy; in that case we build all, i.e. both the constants
  // and the sections
  manager->build(all);

  // We start the VRMLSectionsVisitor to build the vrml file
  // If we don't specify a string in the constructor, it will build the
  // vrml file for all the volumes placed in the topVolume, otherwise it
  // will build the vrml file for the specified volume. The output is
  // placed in sections.vrml
  if (argc == 2)
    manager->startVisitor(new GDDVRMLSectionsVisitor("")); 
  else
    manager->startVisitor(new GDDVRMLSectionsVisitor(argv[2])); 

  // We start the HTMLConstantsVisitor to build the html file with the
  // constants tables. Colors and layout are stolen from Joanne ones.
  manager->startVisitor(new GDDHTMLConstantsVisitor());

  // We retrive the hierarchy entry point, i.e. the GDD object. It
  // contains all the relevant information
  GDD* g = manager->getGDD();
  // We can use some methods of the GDD object to retrive relevant information,
  // for example the value for a constant or the dimension for a volume.
  // In this example we extract the value for the char constant convMat
  // and the material and X dimension for the volume (a box) CsIElementEnv
  std::cout << "The convMat value is " << g->getCharConstByName("convMat") << std::endl;
  std::cout << "The material of the CsIElementEnv is " 
	    << (static_cast<GDDshape*>(g->getVolumeByName("CsIElementEnv")))->getMaterial() << std::endl;
  std::cout << "and its X dimension (in mm) is " 
	    << (static_cast<GDDbox*>(g->getVolumeByName("CsIElementEnv")))->getX() << std::endl;


  delete manager;
  return(0);
}

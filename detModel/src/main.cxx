/// Test program for the generic model

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDVRMLSectionsVisitor.h"
#include "detModel/Management/GDDprinterSectionsVisitor.h"
#include "detModel/Management/GDDHTMLConstantsVisitor.h"
#include "detModel/Management/GDDFakeBuilder.h"
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
  manager->setBuilder(new GDDFakeBuilder);

  manager->setNameFile(argv[1]);
  

  // We set the mode for the choice elements in the XML file
  manager->setMode("digi recon");

  // We build the hierarchy; in that case we build all, i.e. both the constants
  // and the sections
  manager->build(all);

  // We start the VRMLSectionsVisitor to build the vrml file
  // If we don't specify a string in the constructor, it will build the
  // vrml file for all the volumes placed in the topVolume, otherwise it
  // will build the vrml file for the specified volume. The output is
  // placed in sections.vrml
  GDDVRMLSectionsVisitor* visitor;
  if (argc == 2)
    visitor = new GDDVRMLSectionsVisitor("");  
  else
    visitor = new GDDVRMLSectionsVisitor(argv[2]);  

  visitor->setOpacity("FOAM05",0.0);
  // visitor->setDepth("oneTKR", 0);
  // visitor->setDepth("oneCAL", 0);
  manager->startVisitor(visitor);
  
  // We start the HTMLConstantsVisitor to build the html file with the
  // constants tables. Colors and layout are stolen from Joanne ones.
  manager->startVisitor(new GDDHTMLConstantsVisitor());

  // We retrive the hierarchy entry point, i.e. the GDD object. It
  // contains all the relevant information
  GDD* g = manager->getGDD();

  delete manager;
  return(0);
}






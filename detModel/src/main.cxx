/// Test program for the generic model

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDVRMLSectionsVisitor.h"
#include "detModel/Management/GDDprinterSectionsVisitor.h"
#include "detModel/Management/GDDXercesBuilder.h"

int main(int argc, char* argv[]) {
  
  if (argc == 1) return 0;
  
  GDDmanager* manager = GDDmanager::getPointer();
  
  manager->setBuilder(new GDDXercesBuilder(argv[1]));
  manager->setMode("digi recon");
  manager->build();

  if (argc == 2)
    manager->startVisitor(new GDDVRMLSectionsVisitor()); 
  else
    manager->startVisitor(new GDDVRMLSectionsVisitor(argv[2])); 
  // manager->startVisitor(new GDDprinterSectionsVisitor);

  delete manager;
  return(0);
}

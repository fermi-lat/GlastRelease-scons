#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Gdd.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Materials/Element.h"
#include "detModel/Materials/Composite.h"
#include "detModel/Management/PrinterMaterialsVisitor.h"
#include "detModel/Management/Manager.h"

namespace detModel{

  PrinterMaterialsVisitor::PrinterMaterialsVisitor()
  {
    setRecursive(1);
    m_gdd = 0;
  };

  void PrinterMaterialsVisitor::visitGdd(Gdd* gdd)
  {
    std::cout << "Starting the visit of the Gdd" << std::endl;
    std::cout << "DTDversion: " << gdd->getDTDversion() << std::endl;
    std::cout << "CVSid: " << gdd->getCVSid() << std::endl;
    std::cout << "-----------------------------" << std::endl;
    std::cout << " " << std::endl;
    m_gdd = gdd;
  }
  
  void  PrinterMaterialsVisitor::visitMatCollection(MatCollection* materials)
  {
    std::cout << "\t Starting the visit of materials list version " << materials->getVersion() << std::endl;
    std::cout << "\t written by " << materials->getAuthor()  << " (" << materials->getDate() << ")" << std::endl;
    std::cout << "\t composed by " << materials->getMaterials().size() << " materials." << std::endl;
    std::cout << "\t -------------------------------------------------------------------" << std::endl;
    std::cout << " " << std::endl;
  }

  void  PrinterMaterialsVisitor::visitElement(Element* element)
  {
    std::cout << "\t\t Starting the visit of element " << element->getName() << std::endl;  

    if (element->getSymbol() != "")
      std::cout << "\t\t whose symbol is " << element->getSymbol() << std::endl;

    std::cout << "\t\t whose z is " << element->getZ() << std::endl;

    if (element->getDensity())
      std::cout << "\t\t whose density is " << element->getDensity() << std::endl;
    std::cout << "\t\t whose aweight is " << element->getAweight() << std::endl;

    //    Manager* man = Manager::getPointer();
    //    Gdd* g = man->getGdd();

    if(Color* col = m_gdd->getMaterialColorByName(element->getName()))
      {
	std::cout << "\t\t whose color is (RGBT) ";
	std::cout << col->getRed() << " " <<
	  col->getGreen() << " " << col->getBlue() << " " << col->getTra() << std::endl;
      }
    
    std::cout << "\t\t -----------------------------------" << std::endl;
  }


  void  PrinterMaterialsVisitor::visitComposite(Composite* comp)
  {
    unsigned int i;

    std::cout << "\t\t Starting the visit of element " << comp->getName() << std::endl;      
    std::cout << "\t\t whose density is " << comp->getDensity() << std::endl;    

    std::cout << "\t\t composed by " << comp->getNComponents() << " components" << std::endl;

    for(i=0;i<comp->getComponents().size();i++)
      {
	std::cout << "\t\t\t -- "<< comp->getComponents()[i]->getName();
	if (comp->isFractions())
	  std::cout << " -> fract. = " << comp->getFractions()[i] << std::endl;
	else
	  std::cout << " -> n atoms = " << comp->getAtoms()[i] << std::endl;
      }

    //    Manager* man = Manager::getPointer();
    //    Gdd* g = man->getGdd();

    if(Color* col = m_gdd->getMaterialColorByName(comp->getName()))
      {
	std::cout << "\t\t whose color is (RGBT) ";
	std::cout << col->getRed() << " " <<
	  col->getGreen() << " " << col->getBlue() << " " << col->getTra() << std::endl;
      }

    std::cout << "\t\t -----------------------------------" << std::endl;
  }
}

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/GDD.h"
#include "detModel/Constants/GDDconstants.h"
#include "detModel/Constants/GDDconstCategory.h"
#include "detModel/Constants/GDDintConst.h"
#include "detModel/Constants/GDDfloatConst.h"
#include "detModel/Constants/GDDstringConst.h"
#include "detModel/Constants/GDDdoubleConst.h"

#include "detModel/Management/GDDprinterConstantsVisitor.h"

GDDprinterConstantsVisitor::GDDprinterConstantsVisitor()
{
  setRecursive(1);
};

void GDDprinterConstantsVisitor::visitGDD(GDD* gdd)
{
  std::cout << "Starting the visit of the GDD" << std::endl;
  std::cout << "DTDversion: " << gdd->getDTDversion() << std::endl;
  std::cout << "CVSid: " << gdd->getCVSid() << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << " " << std::endl;
}
  
void  GDDprinterConstantsVisitor::visitConstants(GDDconstants* constants)
{
  std::cout << "\t Starting the visit of constants version " << constants->getVersion() << std::endl;
  std::cout << "\t composed by " << constants->getCategories().size() << " categories." << std::endl;
  std::cout << "\t -------------------------------------------------------------------" << std::endl;
  std::cout << " " << std::endl;
}

void  GDDprinterConstantsVisitor::visitConstCategory(GDDconstCategory* category)
{
  std::cout << "\t\t Starting the visit of category " << category->getName() << std::endl;  
  std::cout << "\t\t that is ";
  if(category->getPrimary())
    std::cout << "primary." << std::endl;
  else
    std::cout << "derived." << std::endl;
  std::cout << "\t\t "<< category->getOverview() << std::endl;  
  std::cout << "\t\t -----------------------------------" << std::endl;
}


void  GDDprinterConstantsVisitor::visitIntConst(GDDintConst* c)
{
  std::cout << "\t\t\t Starting the visit of the int const " << c->getName() << std::endl;  
  std::cout << "\t\t\t " << c->getNote() << std::endl; 
  std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
  std::cout << "\t\t\t -----------------------------------" << std::endl;
}

void  GDDprinterConstantsVisitor::visitFloatConst(GDDfloatConst* c)
{
  std::cout << "\t\t\t Starting the visit of the float const " << c->getName() << std::endl;  
  std::cout << "\t\t\t " << c->getNote() << std::endl; 
  std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
  std::cout << "\t\t\t -----------------------------------" << std::endl;
}

void  GDDprinterConstantsVisitor::visitDoubleConst(GDDdoubleConst* c)
{
  std::cout << "\t\t\t Starting the visit of the double const " << c->getName() << std::endl;  
  std::cout << "\t\t\t " << c->getNote() << std::endl; 
  std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
  std::cout << "\t\t\t -----------------------------------" << std::endl;
}

void  GDDprinterConstantsVisitor::visitStringConst(GDDstringConst* c)
{
  std::cout << "\t\t\t Starting the visit of the string const " << c->getName() << std::endl;  
  std::cout << "\t\t\t " << c->getNote() << std::endl; 
  std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
  std::cout << "\t\t\t -----------------------------------" << std::endl;
}













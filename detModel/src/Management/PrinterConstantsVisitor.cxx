#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Gdd.h"
#include "detModel/Constants/Constants.h"
#include "detModel/Constants/ConstCategory.h"
#include "detModel/Constants/IntConst.h"
#include "detModel/Constants/FloatConst.h"
#include "detModel/Constants/StringConst.h"
#include "detModel/Constants/DoubleConst.h"
#include "detModel/Management/PrinterConstantsVisitor.h"

namespace detModel{

  PrinterConstantsVisitor::PrinterConstantsVisitor()
  {
    setRecursive(1);
  };

  void PrinterConstantsVisitor::visitGdd(Gdd* gdd)
  {
    std::cout << "Starting the visit of the Gdd" << std::endl;
    std::cout << "DTDversion: " << gdd->getDTDversion() << std::endl;
    std::cout << "CVSid: " << gdd->getCVSid() << std::endl;
    std::cout << "-----------------------------" << std::endl;
    std::cout << " " << std::endl;
  }
  
  void  PrinterConstantsVisitor::visitConstants(Constants* constants)
  {
    std::cout << "\t Starting the visit of constants version " << constants->getVersion() << std::endl;
    std::cout << "\t composed by " << constants->getCategories().size() << " categories." << std::endl;
    std::cout << "\t -------------------------------------------------------------------" << std::endl;
    std::cout << " " << std::endl;
  }

  void  PrinterConstantsVisitor::visitConstCategory(ConstCategory* category)
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


  void  PrinterConstantsVisitor::visitIntConst(IntConst* c)
  {
    std::cout << "\t\t\t Starting the visit of the int const " << c->getName() << std::endl;  
    std::cout << "\t\t\t " << c->getNote() << std::endl; 
    std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
    std::cout << "\t\t\t -----------------------------------" << std::endl;
  }

  void  PrinterConstantsVisitor::visitFloatConst(FloatConst* c)
  {
    std::cout << "\t\t\t Starting the visit of the float const " << c->getName() << std::endl;  
    std::cout << "\t\t\t " << c->getNote() << std::endl; 
    std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
    std::cout << "\t\t\t -----------------------------------" << std::endl;
  }

  void  PrinterConstantsVisitor::visitDoubleConst(DoubleConst* c)
  {
    std::cout << "\t\t\t Starting the visit of the double const " << c->getName() << std::endl;  
    std::cout << "\t\t\t " << c->getNote() << std::endl; 
    std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
    std::cout << "\t\t\t -----------------------------------" << std::endl;
  }

  void  PrinterConstantsVisitor::visitStringConst(StringConst* c)
  {
    std::cout << "\t\t\t Starting the visit of the string const " << c->getName() << std::endl;  
    std::cout << "\t\t\t " << c->getNote() << std::endl; 
    std::cout << "\t\t\t whose value is " << c->getValue() << std::endl;
    std::cout << "\t\t\t -----------------------------------" << std::endl;
  }

}











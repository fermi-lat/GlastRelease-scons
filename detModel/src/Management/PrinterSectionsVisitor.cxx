#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "detModel/Gdd.h"
#include "detModel/Sections/Section.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/Ensemble.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Sections/IdField.h"
#include "detModel/Sections/Seg.h"
#include "detModel/Management/PrinterSectionsVisitor.h"

namespace detModel{

  PrinterSectionsVisitor::PrinterSectionsVisitor()
  {
    setRecursive(1);
  };

  void PrinterSectionsVisitor::visitGdd(Gdd* gdd)
  {
    std::cout << "Starting the visit of the Gdd" << std::endl;
    std::cout << "DTDversion: " << gdd->getDTDversion() << std::endl;
    std::cout << "CVSid: " << gdd->getCVSid() << std::endl;
    std::cout << "-----------------------------" << std::endl;
    std::cout << " " << std::endl;
  }
  
  void  PrinterSectionsVisitor::visitSection(Section* section)
  {
    std::cout << "\t Starting the visit of section " << section->getName() << std::endl;
    std::cout << "\t composed by " << section->getVolumes().size() << " volumes and " 
	      << section->getChoices().size() << " choices" << std::endl;
    std::cout << "\t whose topVolume is " << section->getTopVolume()->getName() << std::endl;
    std::cout << "\t Author: " << section->getAuthor() << std::endl;
    std::cout << "\t Version: " << section->getVersion() << std::endl;
    std::cout << "\t Date: " << section->getDate() << std::endl;
    std::cout << "\t -----------------------------------" << std::endl;
    std::cout << " " << std::endl;
  }

  void  PrinterSectionsVisitor::visitEnsemble(Ensemble* ensemble)
  {
    if(Composition* comp = dynamic_cast<Composition*>(ensemble))
      {
	std::cout << "\t\t Starting the visit of composition " << comp->getName() << std::endl;  
	std::cout << "\t\t whose envelope is " << comp->getEnvelope()->getName() << std::endl;  
	std::cout << "\t\t -----------------------------------" << std::endl;
      }
    else if(Stack* st = dynamic_cast<Stack*>(ensemble))
      {
	std::cout << "\t\t Starting the visit of stack " << st->getName() << std::endl;  
	std::cout << "\t\t that is oriented along ";
	switch(st->getAxisDir()){
	case (Stack::xDir):
	  std::cout << "the X axis";
	  break;
	case (Stack::yDir):
	  std::cout << "the Y axis";
	  break;
	case (Stack::zDir):
	  std::cout << "the Z axis";
	  break;    
	}   
      }

    std::cout << " " << std::endl;
    std::cout << "\t\t -----------------------------------" << std::endl;
  }

  void  PrinterSectionsVisitor::visitBox(Box* box)
  {
    std::cout << "\t\t Visiting the box " << box->getName() << std::endl;  
    std::cout << "\t\t whose material is " << box->getMaterial() << std::endl;  
    std::cout << "\t\t whose dimensions are   X: " 
	      << box->getX() << "  Y: " << box->getY() << "  Z: " << box->getZ() << std::endl;  
    std::cout << "\t\t --------------------------------------------------------" << std::endl;
  }

  void  PrinterSectionsVisitor::visitPosXYZ(PosXYZ* pos)
  {
    std::cout << "\t\t\t Visiting the PosXYZ of " << pos->getVolume()->getName() << std::endl;  
    std::cout << "\t\t\t whose position is   X: " 
	      << pos->getX() << "  Y: " << pos->getY() << "  Z: " << pos->getZ() << std::endl;  
    std::cout << "\t\t\t ------------------ -----------------------------------" << std::endl;
  }

  void  PrinterSectionsVisitor::visitAxisMPos(AxisMPos* pos)
  {
    std::cout << "\t\t\t Visiting the AxisMPos of " << pos->getVolume()->getName() << std::endl;    
    std::cout << "\t\t\t that is replicated " << pos->getNcopy() << " times" << std::endl;
    std::cout << "\t\t\t with displacements along the axis by " << std::endl;
    for(unsigned int i = 0; i<pos->getNcopy();i++)
      std::cout << "\t\t\t\t  " << pos->getDisp(i) << std::endl;
    std::cout << "\t\t\t -----------------------------------" << std::endl;
  }

  void  PrinterSectionsVisitor::visitIdField(IdField* idf)
  {
    std::cout << "\t\t\t\t Visiting the idField " << idf->getName() << std::endl;    
    std::cout << "\t\t\t\t whose value is " << idf->getValue() << " and step is " <<
      idf->getStep() << std::endl;
    std::cout << "\t\t\t\t -----------------------------------" << std::endl;
  }



  void  PrinterSectionsVisitor::visitSeg(Seg* s)
  {
    std::cout << "\t\t\t\t Visiting seg oriented in direction  : " << (int)(s->getAxis()) << std::endl;    
    std::cout << "\t\t\t\t whose reason is " << s->getReason() << "  nSeg is " <<
      s->getnSeg()  <<std::endl;
    std::cout << "\t\t\t\t -----------------------------------" << std::endl;
  }

}














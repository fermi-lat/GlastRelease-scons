#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "detModel/GDD.h"
#include "detModel/Sections/GDDsection.h"
#include "detModel/Sections/GDDbox.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDensemble.h"
#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDaxisMPos.h"
#include "detModel/Sections/GDDidField.h"
#include "detModel/Sections/GDDseg.h"
#include "detModel/Management/GDDprinterSectionsVisitor.h"

GDDprinterSectionsVisitor::GDDprinterSectionsVisitor()
{
  setRecursive(1);
};

void GDDprinterSectionsVisitor::visitGDD(GDD* gdd)
{
  std::cout << "Starting the visit of the GDD" << std::endl;
  std::cout << "DTDversion: " << gdd->getDTDversion() << std::endl;
  std::cout << "CVSid: " << gdd->getCVSid() << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << " " << std::endl;
}
  
void  GDDprinterSectionsVisitor::visitSection(GDDsection* section)
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

void  GDDprinterSectionsVisitor::visitEnsemble(GDDensemble* ensemble)
{
  if(GDDcomposition* comp = dynamic_cast<GDDcomposition*>(ensemble))
    {
      std::cout << "\t\t Starting the visit of composition " << comp->getName() << std::endl;  
      std::cout << "\t\t whose envelope is " << comp->getEnvelope()->getName() << std::endl;  
      std::cout << "\t\t -----------------------------------" << std::endl;
    }
  else if(GDDstack* st = dynamic_cast<GDDstack*>(ensemble))
    {
      std::cout << "\t\t Starting the visit of stack " << st->getName() << std::endl;  
      std::cout << "\t\t that is oriented along ";
      switch(st->getAxisDir()){
      case (GDDstack::xDir):
	std::cout << "the X axis";
	break;
      case (GDDstack::yDir):
	std::cout << "the Y axis";
	break;
      case (GDDstack::zDir):
	std::cout << "the Z axis";
	break;    
      }   
    }

  std::cout << " " << std::endl;
  std::cout << "\t\t -----------------------------------" << std::endl;
}

void  GDDprinterSectionsVisitor::visitBox(GDDbox* box)
{
  std::cout << "\t\t Visiting the box " << box->getName() << std::endl;  
  std::cout << "\t\t whose material is " << box->getMaterial() << std::endl;  
  std::cout << "\t\t whose dimensions are   X: " 
	    << box->getX() << "  Y: " << box->getY() << "  Z: " << box->getZ() << std::endl;  
  std::cout << "\t\t --------------------------------------------------------" << std::endl;
}

void  GDDprinterSectionsVisitor::visitPosXYZ(GDDposXYZ* pos)
{
  std::cout << "\t\t\t Visiting the PosXYZ of " << pos->getVolume()->getName() << std::endl;  
  std::cout << "\t\t\t whose position is   X: " 
	    << pos->getX() << "  Y: " << pos->getY() << "  Z: " << pos->getZ() << std::endl;  
  std::cout << "\t\t\t ------------------ -----------------------------------" << std::endl;
}

void  GDDprinterSectionsVisitor::visitAxisMPos(GDDaxisMPos* pos)
{
  std::cout << "\t\t\t Visiting the AxisMPos of " << pos->getVolume()->getName() << std::endl;    
  std::cout << "\t\t\t that is replicated " << pos->getNcopy() << " times" << std::endl;
  std::cout << "\t\t\t with displacements along the axis by " << std::endl;
  for(unsigned int i = 0; i<pos->getNcopy();i++)
    std::cout << "\t\t\t\t  " << pos->getDisp(i) << std::endl;
  std::cout << "\t\t\t -----------------------------------" << std::endl;
}

void  GDDprinterSectionsVisitor::visitIdField(GDDidField* idf)
{
  std::cout << "\t\t\t\t Visiting the idField " << idf->getName() << std::endl;    
  std::cout << "\t\t\t\t whose value is " << idf->getValue() << " and step is " <<
    idf->getStep() << std::endl;
  std::cout << "\t\t\t\t -----------------------------------" << std::endl;
}



void  GDDprinterSectionsVisitor::visitSeg(GDDseg* s)
{
  std::cout << "\t\t\t\t Visiting seg oriented in direction  : " << (int)(s->getAxis()) << std::endl;    
  std::cout << "\t\t\t\t whose reason is " << s->getReason() << "  nSeg is " <<
    s->getnSeg()  <<std::endl;
  std::cout << "\t\t\t\t -----------------------------------" << std::endl;
}
















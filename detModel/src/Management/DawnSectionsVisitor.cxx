#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Materials/MatCollection.h"
#include "detModel/Management/DawnSectionsVisitor.h"
#include "detModel/Management/Manager.h"
#include "detModel/Sections/Section.h"
#include "detModel/Gdd.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Sections/IdField.h"
#include "detModel/Sections/Position.h"
#include "detModel/Utilities/ColorCreator.h"
#include "detModel/Utilities/Color.h"
#include "detModel/Utilities/Vector.h"

namespace detModel{

DawnSectionsVisitor::DawnSectionsVisitor(std::string nvol)
{
  //  typedef map<std::string,Color*>M1;

  setRecursive(0);
  actualVolume = nvol;

  out.open("sections.prim");

  Manager* manager = Manager::getPointer();
  Gdd* g = manager->getGdd();

  colorsMap = g->getMaterials()->getMaterialColors();
  
  actualPos = new Vector(0,0,0);
  actualXAxis = new Vector(1,0,0);
  actualYAxis = new Vector(0,1,0);
  actualZAxis = new Vector(0,0,1);

  actualXRot = 0;
  actualYRot = 0;
  actualZRot = 0;

};

DawnSectionsVisitor::~DawnSectionsVisitor()
{
  out.close();
}


void DawnSectionsVisitor::visitGdd(Gdd* Gdd)
{
  typedef std::vector<Section*>sec;
  std::vector <Section*>::iterator i;

  sec s = Gdd->getSections();
  for(i=s.begin(); i!=s.end();i++)
    (*i)->AcceptNotRec(this);
}
  
void  DawnSectionsVisitor::visitSection(Section* section)
{
  Volume* vol=0;
  std::vector <Volume*>::iterator v;
  
  if (actualVolume == "")
    {
      std::vector<Volume*> volumes = section->getVolumes();
      
      for(v=volumes.begin(); v!=volumes.end(); v++){
	if(Ensemble* ens = dynamic_cast<Ensemble*>(*v))
	  {
	    if (ens == section->getTopVolume())
	      vol = ens;
	  }
      }
    }
  else
    {
      Manager* manager = Manager::getPointer();
      if (manager->getGdd()->getVolumeByName(actualVolume))
	vol = manager->getGdd()->getVolumeByName(actualVolume);
      else
	{
	  std::cout << "No such volume" << std::endl;
	  exit(0);
	}
    }

  /// Here we start the visit of the hierarchy.
  out << "##G4.PRIM-FORMAT-2.4" << std::endl;
  out << "/BoundingBox ";
  out << - vol->getBBox()->getXDim()/2 << " ";
  out << - vol->getBBox()->getYDim()/2 << " ";
  out << - vol->getBBox()->getZDim()/2 << " ";
  out << vol->getBBox()->getXDim()/2 << " ";
  out << vol->getBBox()->getYDim()/2 << " ";
  out << vol->getBBox()->getZDim()/2 << std::endl;

  out << "!SetCamera" << std::endl;
  out << "!OpenDevice" << std::endl;
  out << "!BeginModeling" << std::endl;          

  out << "/ForceWireframe 1" << std::endl;
  out << "/Box " <<
    vol->getBBox()->getXDim()/2 << " " <<
    vol->getBBox()->getYDim()/2 << " " <<
    vol->getBBox()->getZDim()/2 << std::endl;
  out << "/ForceWireframe 0" << std::endl;

  out << "/ColorRGB 0 0 0" << std::endl; 
  out << "/FontName Times-Roman" << std::endl; 
  out << "/Text2DS 30 40 12 " << vol->getName() << "(";
  out << vol->getBBox()->getXDim() << " ";
  out << vol->getBBox()->getYDim() << " ";
  out << vol->getBBox()->getZDim() << ")" << std::endl;


  vol->AcceptNotRec(this);
  out << "!EndModeling" << std::endl;
  out << "!DrawAll" << std::endl;
  out << "!CloseDevice" << std::endl;          
}

void  DawnSectionsVisitor::visitEnsemble(Ensemble* ensemble)
{
  std::vector <Position*>::iterator i;
  typedef std::vector<Position*> pos;
  
  out << "# " << ensemble->getName() << std::endl;
      
  /// Here the positioned volumes are visited
  pos p = ensemble->getPositions();
  for(i=p.begin(); i!=p.end();i++)
    (*i)->AcceptNotRec(this);

      
  /// Here the envelope is visited if the ensamble is a composition
  if (Composition* comp = dynamic_cast<Composition*>(ensemble))
    {
      out << "/ForceWireframe 1" << std::endl;
      comp->getEnvelope()->AcceptNotRec(this);
      out << "/ForceWireframe 0" << std::endl;
    }
}


void  DawnSectionsVisitor::visitBox(Box* box)
{
  typedef std::map<std::string, Color*> M;
  M::const_iterator j; 

  j = colorsMap.find(box->getMaterial());
  if (j == colorsMap.end()) return;


  if(j->second->getTra() != 0)
    out << "/ForceWireframe 1" << std::endl;

  out << "/Origin " << actualPos->getX() << " " << 
    actualPos->getY() << " " << actualPos->getZ() << std::endl;

  out << "/BaseVector " << 
    actualXAxis->getX() << " " << 
    actualXAxis->getY() << " " << 
    actualXAxis->getZ() << " " <<
    actualYAxis->getX() << " " << 
    actualYAxis->getY() << " " << 
    actualYAxis->getZ() << std::endl;

  
  out << "/ColorRGB " << j->second->getRed() << " " 
      <<  j->second->getGreen() << " " 
      <<  j->second->getBlue()  << std::endl;
    

  out << "/Box " <<
    box->getX()/2 << " " << box->getY()/2 << " " << box->getZ()/2 << " " << std::endl;

  out << "/ForceWireframe 0" << std::endl;
}

void  DawnSectionsVisitor::visitPosXYZ(PosXYZ* pos)
{
  Vector* temp = new Vector(actualPos->getX(),
				  actualPos->getY(),
				  actualPos->getZ()); 

  Vector* tempX = new Vector(actualXAxis->getX(),
				   actualXAxis->getY(),
				   actualXAxis->getZ()); 

  Vector* tempY = new Vector(actualYAxis->getX(),
				   actualYAxis->getY(),
				   actualYAxis->getZ()); 

  Vector* tempZ = new Vector(actualZAxis->getX(),
				   actualZAxis->getY(),
				   actualZAxis->getZ()); 


  Vector* relPos = new Vector(pos->getX(),pos->getY(),pos->getZ());

  actualXAxis->Rotate(0,pos->getXRot()); 
  actualYAxis->Rotate(0,pos->getXRot()); 
  actualZAxis->Rotate(0,pos->getXRot()); 
  
  actualXAxis->Rotate(1,pos->getYRot()); 
  actualYAxis->Rotate(1,pos->getYRot()); 
  actualZAxis->Rotate(1,pos->getYRot()); 

  actualXAxis->Rotate(2,pos->getZRot()); 
  actualYAxis->Rotate(2,pos->getZRot()); 
  actualZAxis->Rotate(2,pos->getZRot()); 

  
  relPos->Rotate(0,actualXRot);
  relPos->Rotate(1,actualYRot);
  relPos->Rotate(2,actualZRot);

  
  actualPos->x = temp->x + relPos->x;
  actualPos->y = temp->y + relPos->y;
  actualPos->z = temp->z + relPos->z;

  pos->getVolume()->AcceptNotRec(this);

  actualPos->x = temp->x;
  actualPos->y = temp->y;
  actualPos->z = temp->z;

  actualXAxis->x = tempX->x;
  actualXAxis->y = tempX->y;
  actualXAxis->z = tempX->z;

  actualYAxis->x = tempY->x;
  actualYAxis->y = tempY->y;
  actualYAxis->z = tempY->z;

  actualZAxis->x = tempZ->x;
  actualZAxis->y = tempZ->y;
  actualZAxis->z = tempZ->z;

  delete temp;
  delete tempX;
  delete tempY;
  delete tempZ;
}


void  DawnSectionsVisitor::visitAxisMPos(AxisMPos* pos)
{
  int i;
  int n;

  double tempXRot = actualXRot;
  double tempYRot = actualYRot;
  double tempZRot = actualZRot;

  Vector* temp = new Vector(0,0,0);

  Vector* relPos = new Vector(0,0,0);

  Vector* tempX = new Vector(actualXAxis->getX(),
				   actualXAxis->getY(),
				   actualXAxis->getZ()); 

  Vector* tempY = new Vector(actualYAxis->getX(),
				   actualYAxis->getY(),
				   actualYAxis->getZ()); 

  Vector* tempZ = new Vector(actualZAxis->getX(),
				   actualZAxis->getY(),
				   actualZAxis->getZ()); 


  temp->x = actualPos->getX();
  temp->y = actualPos->getY();
  temp->z = actualPos->getZ(); 


  n = pos->getNcopy();

  
  for(i=0;i<n;i++){
  switch(pos->getAxisDir()){
  case (Stack::xDir):
    {
      relPos->x = pos->getDx()+ pos->getDisp(i);
      relPos->y = pos->getDy();
      relPos->z = pos->getDz();      

      actualXRot = actualXRot + pos->getRotation();

      actualXAxis->Rotate(0,pos->getRotation()); 
      actualYAxis->Rotate(0,pos->getRotation()); 
      actualZAxis->Rotate(0,pos->getRotation());
 
    }
    break;
  case (Stack::yDir):
    {
      relPos->y = pos->getDy()+ pos->getDisp(i);
      relPos->x = pos->getDx();
      relPos->z = pos->getDz();      

      actualYRot = actualYRot + pos->getRotation();

      actualXAxis->Rotate(1,pos->getRotation()); 
      actualYAxis->Rotate(1,pos->getRotation()); 
      actualZAxis->Rotate(1,pos->getRotation());

    }
    break;
  case (Stack::zDir):
    {
      relPos->z = pos->getDz()+ pos->getDisp(i);
      relPos->y = pos->getDy();
      relPos->x = pos->getDx();      

      actualZRot = actualZRot + pos->getRotation();

      actualXAxis->Rotate(2,pos->getRotation()); 
      actualYAxis->Rotate(2,pos->getRotation()); 
      actualZAxis->Rotate(2,pos->getRotation());

    }
    break;
  }

  relPos->Rotate(0,actualXRot);
  relPos->Rotate(1,actualYRot);
  relPos->Rotate(2,actualZRot);


  actualPos->x = actualPos->x + relPos->x;
  actualPos->y = actualPos->y + relPos->y;
  actualPos->z = actualPos->z + relPos->z;

  pos->getVolume()->AcceptNotRec(this);

  actualPos->x = temp->x;
  actualPos->y = temp->y;
  actualPos->z = temp->z;

  actualXAxis->x = tempX->x;
  actualXAxis->y = tempX->y;
  actualXAxis->z = tempX->z;

  actualYAxis->x = tempY->x;
  actualYAxis->y = tempY->y;
  actualYAxis->z = tempY->z;

  actualZAxis->x = tempZ->x;
  actualZAxis->y = tempZ->y;
  actualZAxis->z = tempZ->z;

  actualXRot = tempXRot;
  actualYRot = tempYRot;
  actualZRot = tempZRot;
  }

  delete temp;
  delete tempX;
  delete tempY;
  delete tempZ;
}

void  DawnSectionsVisitor::visitIdField(IdField*)
{

}

void  DawnSectionsVisitor::visitSeg(Seg*)
{

}



}

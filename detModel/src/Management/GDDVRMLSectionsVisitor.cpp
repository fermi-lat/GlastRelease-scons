#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Management/GDDVRMLSectionsVisitor.h"
#include "detModel/Management/GDDmanager.h"
#include "detModel/Sections/GDDsection.h"
#include "detModel/GDD.h"
#include "detModel/Sections/GDDbox.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDaxisPos.h"
#include "detModel/Sections/GDDaxisMPos.h"
#include "detModel/Sections/GDDidField.h"
#include "detModel/Sections/GDDanyPosition.h"

GDDVRMLSectionsVisitor::GDDVRMLSectionsVisitor()
{
  setRecursive(0);
  setType(sectionsVisitor);
  actualVolume = "";
};

GDDVRMLSectionsVisitor::GDDVRMLSectionsVisitor(string nvol)
{
  setRecursive(0);
  setType(sectionsVisitor);
  actualVolume = nvol;
};

void GDDVRMLSectionsVisitor::visitGDD(GDD* gdd)
{
  unsigned int i;
  std::cout << "#VRML V2.0 utf8 " << std::endl;

  std::cout << " DEF SI " << std::endl;
  std::cout << " Appearance { " <<  std::endl;
  std::cout << " material Material { " <<  std::endl;
  std::cout << "       ambientIntensity  0.5" <<  std::endl;        
  std::cout << "       diffuseColor      0.8 0.8 0.2" <<  std::endl;
  std::cout << "       emissiveColor     0 0 0"  <<  std::endl;     
  std::cout << "       shininess         0.2" <<  std::endl;        
  std::cout << "       specularColor     0 0 0.2" <<  std::endl;      
  std::cout << "       transparency      0.0" <<  std::endl;          
  std::cout << "     }" <<  std::endl;
  std::cout << " }" <<  std::endl;

  std::cout << " DEF PB " << std::endl;
  std::cout << " Appearance { " <<  std::endl;
  std::cout << " material Material { " <<  std::endl;
  std::cout << "       ambientIntensity  0.5" <<  std::endl;        
  std::cout << "       diffuseColor      0.5 0.5 0.5" <<  std::endl;
  std::cout << "       emissiveColor     0 0 0"  <<  std::endl;     
  std::cout << "       shininess         0.2" <<  std::endl;        
  std::cout << "       specularColor     0 0 0.2" <<  std::endl;      
  std::cout << "       transparency      0.3" <<  std::endl;          
  std::cout << "     }" <<  std::endl;
  std::cout << " }" <<  std::endl;

  std::cout << " DEF C " << std::endl;
  std::cout << " Appearance { " <<  std::endl;
  std::cout << " material Material { " <<  std::endl;
  std::cout << "       ambientIntensity  0.5" <<  std::endl;        
  std::cout << "       diffuseColor      0.1 0.1 0.1" <<  std::endl;
  std::cout << "       emissiveColor     0 0 0.2"  <<  std::endl;     
  std::cout << "       shininess         0.2" <<  std::endl;        
  std::cout << "       specularColor     0 0 0.2" <<  std::endl;      
  std::cout << "       transparency      0.8" <<  std::endl;          
  std::cout << "     }" <<  std::endl;
  std::cout << " }" <<  std::endl;

  std::cout << " DEF RUBBER60A1 " << std::endl;
  std::cout << " Appearance { " <<  std::endl;
  std::cout << " material Material { " <<  std::endl;
  std::cout << "       ambientIntensity  0.5" <<  std::endl;        
  std::cout << "       diffuseColor      0.8 0.2 0.2" <<  std::endl;
  std::cout << "       emissiveColor     0 0 0"  <<  std::endl;     
  std::cout << "       shininess         0.2" <<  std::endl;        
  std::cout << "       specularColor     0 0 0" <<  std::endl;      
  std::cout << "       transparency      0.3" <<  std::endl;          
  std::cout << "     }" <<  std::endl;
  std::cout << " }" <<  std::endl;

  std::cout << " DEF FOAM05 " << std::endl;
  std::cout << " Appearance { " <<  std::endl;
  std::cout << " material Material { " <<  std::endl;
  std::cout << "       ambientIntensity  0.5" <<  std::endl;        
  std::cout << "       diffuseColor      0.8 0.8 0.8" <<  std::endl;
  std::cout << "       emissiveColor     0 0 0"  <<  std::endl;     
  std::cout << "       shininess         0.2" <<  std::endl;        
  std::cout << "       specularColor     0 0 0" <<  std::endl;      
  std::cout << "       transparency      0.3" <<  std::endl;          
  std::cout << "     }" <<  std::endl;
  std::cout << " }" <<  std::endl;

  std::cout << " DEF CSI " << std::endl;
  std::cout << " Appearance { " <<  std::endl;
  std::cout << " material Material { " <<  std::endl;
  std::cout << "       ambientIntensity  0.5" <<  std::endl;        
  std::cout << "       diffuseColor      0.2 0.2 0.8" <<  std::endl;
  std::cout << "       emissiveColor     0 0 0"  <<  std::endl;     
  std::cout << "       shininess         0.2" <<  std::endl;        
  std::cout << "       specularColor     0 0 0.2" <<  std::endl;      
  std::cout << "       transparency      0.0" <<  std::endl;          
  std::cout << "     }" <<  std::endl;
  std::cout << " }" <<  std::endl;

  std::cout << " DEF DUMMY " << std::endl;
  std::cout << " Appearance { " <<  std::endl;
  std::cout << " material Material { " <<  std::endl;
  std::cout << "       ambientIntensity  0.5" <<  std::endl;        
  std::cout << "       diffuseColor      1.0 0.0 1.0" <<  std::endl;
  std::cout << "       emissiveColor     0 0 0"  <<  std::endl;     
  std::cout << "       shininess         0.2" <<  std::endl;        
  std::cout << "       specularColor     0 0 0" <<  std::endl;      
  std::cout << "       transparency      0.0" <<  std::endl;          
  std::cout << "     }" <<  std::endl;
  std::cout << " }" <<  std::endl;

  for(i=0; i<gdd->getSections()->size();i++){
    ((* gdd->getSections() )[i])->AcceptNotRec(this);
  }
}
  
void  GDDVRMLSectionsVisitor::visitSection(GDDsection* section)
{
  unsigned int i;

  if (actualVolume == "")
    {
      for(i=0; i<section->getVolumes().size();i++){
	if(section->getVolumes()[i]->getVolumeType() == composition)
	  {
	    GDDcomposition* comp = static_cast<GDDcomposition*>(section->getVolumes()[i]);
	    if (comp->getEnvelope() == section->getTopVolume())
	      {
		section->getVolumes()[i]->AcceptNotRec(this);
	      }
	  }
      }
    }
  else
    {
      GDDmanager* manager = GDDmanager::getPointer();
      
      if (manager->getGDD()->getVolumeByName(actualVolume))
	{
	  manager->getGDD()->getVolumeByName(actualVolume)->AcceptNotRec(this);
	}
    }
}


void  GDDVRMLSectionsVisitor::visitComposition(GDDcomposition* composition)
{
  unsigned int i;

  std::cout << "# " << composition->getName() << std::endl;

  for(i=0; i<composition->getPositions().size();i++){
    GDDanyPosition* pos = composition->getPositions()[i];
    pos->AcceptNotRec(this);
  }
}

void  GDDVRMLSectionsVisitor::visitStack(GDDstack* st)
{
  unsigned int j;
  double deltap, delta;
  GDDanyRelativePosition* apos;

  deltap = 0;

  std::cout << "Transform {  # " << st->getName() << std::endl;
  switch(st->getStackType()){
  case sx:
    std::cout << "translation " << -(st->getBBX())*0.5 << " 0 " << " 0 " << std::endl;
    break;	
  case sy:
    std::cout << "translation " << " 0 " << -(st->getBBY())*0.5 << " 0 " << std::endl;
    break;	
  case sz:
    std::cout << "translation " << " 0 " << " 0 " << -(st->getBBZ())*0.5 << std::endl;
    break;	
  }
  std::cout << "children [ " << std::endl;
  for(j=0;j<st->getPositions().size();j++)
    {
      apos = st->getPositions()[j];
      std::cout << "Transform {  #" << apos->getVolume()->getName() << std::endl;
      switch(st->getStackType()){
      case sx:
	delta = apos->getBBX()/2 + deltap;
	std::cout << "translation " << delta << " 0 " << " 0 " << std::endl;
	deltap += apos->getBBX();
	break;	
      case sy:
	delta = apos->getBBY()/2 + deltap;
	std::cout << "translation " << " 0 " << delta << " 0 " << std::endl;
	deltap += apos->getBBY();
	break;	
      case sz:
	delta = apos->getBBZ()/2 + deltap;
	std::cout << "translation " << " 0 " << " 0 " << delta << std::endl;
	deltap += apos->getBBZ();
	break;	
      }
      std::cout << "children [ " << std::endl;
      apos->AcceptNotRec(this);
      std::cout << "] " << std::endl; 
      std::cout << "} " << std::endl;
    }  

  std::cout << "] " << std::endl; 
  std::cout << "} " << std::endl;
}



void  GDDVRMLSectionsVisitor::visitBox(GDDbox* box)
{
  std::cout << "Shape {   #" << box->getName() << std::endl;

  string mat = box->getMaterial();

  if (mat == "CsI")
    {
      std::cout << "appearance USE CSI" << std::endl;
    }
  else if (mat == "Si")
    {
      std::cout << "appearance USE SI " << std::endl;
    }
  else if (mat == "rubber60A1")
    {
      std::cout << "appearance USE RUBBER60A1 " << std::endl;
    }
  else if (mat == "FOAM05")
    {
      std::cout << "appearance USE FOAM05 " << std::endl;
    }
  else if (mat == "Pb")
    {
      std::cout << "appearance USE PB " << std::endl;
    }
  else if (mat == "C")
    {
      std::cout << "appearance USE C " << std::endl;
    }
  else 
    {
      std::cout << "appearance USE DUMMY " << std::endl;
    }
  std::cout << "  geometry Box { " << std::endl;
  std::cout << "                     size " 
	    << box->getX() << " " << box->getY() << " " << box->getZ() << std::endl; 
  std::cout << "                    }" << std::endl;  
  std::cout << "   }" << std::endl;
}

void  GDDVRMLSectionsVisitor::visitPosXYZ(GDDposXYZ* pos)
{
  std::cout << "Transform { " << std::endl;
  std::cout << "translation " << pos->getX() << " " <<
    pos->getY() << " " << pos->getZ() << 
    std::endl;
  std::cout << "children [ " << std::endl;
  pos->getVolume()->AcceptNotRec(this);
  std::cout << "] " << std::endl; 
  std::cout << "} " << std::endl;
}


void  GDDVRMLSectionsVisitor::visitAxisPos(GDDaxisPos* pos)
{
  std::cout << "Transform { " << std::endl;
  switch(pos->getPosDir()){
  case sx:
    std::cout << "rotation " << " 1 0 0 " <<  pos->getRotation()*3.141927/180 << std::endl;  
    break;
  case sy:
    std::cout << "rotation " << " 0 1 0 " <<  pos->getRotation()*3.141927/180 << std::endl;  
    break;
  case sz:
    std::cout << "rotation " << " 0 0 1 " <<  pos->getRotation()*3.141927/180  << std::endl;  
    break;
  }
  std::cout << "children [ " << std::endl;
  pos->getVolume()->AcceptNotRec(this);
  std::cout << "] " << std::endl; 
  std::cout << "} " << std::endl;

}

void  GDDVRMLSectionsVisitor::visitAxisMPos(GDDaxisMPos* pos)
{
  int i;
  int n;
  float g, w;


  n = pos->getNcopy();
  g = pos->getGap();
	  
  for(i=0;i<n;i++){
    std::cout << "Transform { " << std::endl;
    switch(pos->getPosDir()){
    case sx:
      w = pos->getVolume()->getBBX();	
      std::cout << "translation " << i*(w+g)-(n-1)*(w+g)*0.5  << " 0 " << " 0 " << std::endl;  
      break;
    case sy:
      w = pos->getVolume()->getBBY();
      std::cout << "translation " << " 0 " << i*(w+g)-(n-1)*(w+g)*0.5  << " 0 " << std::endl;  
      break;
    case sz:
      w = pos->getVolume()->getBBZ();
      std::cout << "translation " << " 0 " << " 0 " << i*(w+g)-(n-1)*(w+g)*0.5  << std::endl;  
      break;
    }
    std::cout << "children [ " << std::endl;
    pos->getVolume()->AcceptNotRec(this);
    std::cout << "] " << std::endl; 
    std::cout << "} " << std::endl;
  }


}

void  GDDVRMLSectionsVisitor::visitIdField(GDDidField*)
{

}

void  GDDVRMLSectionsVisitor::visitSeg(GDDseg*)
{

}







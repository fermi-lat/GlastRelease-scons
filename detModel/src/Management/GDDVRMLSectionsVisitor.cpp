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

GDDVRMLSectionsVisitor::GDDVRMLSectionsVisitor(string nvol)
{
  unsigned int i;
  typedef map<string,float>M1;
  typedef map<string,int>M2;
  typedef map<string,GDDvolume*>M3;
  M3::const_iterator j;

  setRecursive(0);
  setType(sectionsVisitor);
  actualVolume = nvol;

  out.open("sections.wrl");

  GDDmanager* manager = GDDmanager::getPointer();
  GDD* g = manager->getGDD();
  
  /// We initialize the opacity map
  for(i=0;i<g->getMaterialNames().size();i++)
    opacityMap.insert(M1::value_type(g->getMaterialNames()[i]->getName(),1.0));  
  /// We initialize the depth map
  for(j=g->getVolumesMap().begin();j!=g->getVolumesMap().end();j++)
    depthMap.insert(M1::value_type(j.first(),200));  
};

GDDVRMLSectionsVisitor::~GDDVRMLSectionsVisitor()
{
  out.close();
}

void GDDVRMLSectionsVisitor::visitGDD(GDD* gdd)
{
  unsigned int i;
  out << "#VRML V2.0 utf8 " << std::endl;

  out << " DEF SI " << std::endl;
  out << " Appearance { " <<  std::endl;
  out << " material Material { " <<  std::endl;
  out << "       ambientIntensity  0.5" <<  std::endl;        
  out << "       diffuseColor      0.8 0.8 0.2" <<  std::endl;
  out << "       emissiveColor     0 0 0"  <<  std::endl;     
  out << "       shininess         0.2" <<  std::endl;        
  out << "       specularColor     0 0 0.2" <<  std::endl;      
  out << "       transparency      0.0" <<  std::endl;          
  out << "     }" <<  std::endl;
  out << " }" <<  std::endl;

  out << " DEF PB " << std::endl;
  out << " Appearance { " <<  std::endl;
  out << " material Material { " <<  std::endl;
  out << "       ambientIntensity  0.5" <<  std::endl;        
  out << "       diffuseColor      0.5 0.5 0.5" <<  std::endl;
  out << "       emissiveColor     0 0 0"  <<  std::endl;     
  out << "       shininess         0.2" <<  std::endl;        
  out << "       specularColor     0 0 0.2" <<  std::endl;      
  out << "       transparency      0.3" <<  std::endl;          
  out << "     }" <<  std::endl;
  out << " }" <<  std::endl;

  out << " DEF C " << std::endl;
  out << " Appearance { " <<  std::endl;
  out << " material Material { " <<  std::endl;
  out << "       ambientIntensity  0.5" <<  std::endl;        
  out << "       diffuseColor      0.1 0.1 0.1" <<  std::endl;
  out << "       emissiveColor     0 0 0.2"  <<  std::endl;     
  out << "       shininess         0.2" <<  std::endl;        
  out << "       specularColor     0 0 0.2" <<  std::endl;      
  out << "       transparency      0.8" <<  std::endl;          
  out << "     }" <<  std::endl;
  out << " }" <<  std::endl;

  out << " DEF RUBBER60A1 " << std::endl;
  out << " Appearance { " <<  std::endl;
  out << " material Material { " <<  std::endl;
  out << "       ambientIntensity  0.5" <<  std::endl;        
  out << "       diffuseColor      0.8 0.2 0.2" <<  std::endl;
  out << "       emissiveColor     0 0 0"  <<  std::endl;     
  out << "       shininess         0.2" <<  std::endl;        
  out << "       specularColor     0 0 0" <<  std::endl;      
  out << "       transparency      0.3" <<  std::endl;          
  out << "     }" <<  std::endl;
  out << " }" <<  std::endl;

  out << " DEF FOAM05 " << std::endl;
  out << " Appearance { " <<  std::endl;
  out << " material Material { " <<  std::endl;
  out << "       ambientIntensity  0.5" <<  std::endl;        
  out << "       diffuseColor      0.8 0.8 0.8" <<  std::endl;
  out << "       emissiveColor     0 0 0"  <<  std::endl;     
  out << "       shininess         0.2" <<  std::endl;        
  out << "       specularColor     0 0 0" <<  std::endl;      
  out << "       transparency      0.3" <<  std::endl;          
  out << "     }" <<  std::endl;
  out << " }" <<  std::endl;

  out << " DEF CSI " << std::endl;
  out << " Appearance { " <<  std::endl;
  out << " material Material { " <<  std::endl;
  out << "       ambientIntensity  0.5" <<  std::endl;        
  out << "       diffuseColor      0.2 0.2 0.8" <<  std::endl;
  out << "       emissiveColor     0 0 0"  <<  std::endl;     
  out << "       shininess         0.2" <<  std::endl;        
  out << "       specularColor     0 0 0.2" <<  std::endl;      
  out << "       transparency      0.0" <<  std::endl;          
  out << "     }" <<  std::endl;
  out << " }" <<  std::endl;

  out << " DEF DUMMY " << std::endl;
  out << " Appearance { " <<  std::endl;
  out << " material Material { " <<  std::endl;
  out << "       ambientIntensity  0.5" <<  std::endl;        
  out << "       diffuseColor      1.0 0.0 1.0" <<  std::endl;
  out << "       emissiveColor     0 0 0"  <<  std::endl;     
  out << "       shininess         0.2" <<  std::endl;        
  out << "       specularColor     0 0 0" <<  std::endl;      
  out << "       transparency      0.0" <<  std::endl;          
  out << "     }" <<  std::endl;
  out << " }" <<  std::endl;

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

  out << "# " << composition->getName() << std::endl;

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

  out << "Transform {  # " << st->getName() << std::endl;
  switch(st->getStackType()){
  case sx:
    out << "translation " << -(st->getBBX())*0.5 << " 0 " << " 0 " << std::endl;
    break;	
  case sy:
    out << "translation " << " 0 " << -(st->getBBY())*0.5 << " 0 " << std::endl;
    break;	
  case sz:
    out << "translation " << " 0 " << " 0 " << -(st->getBBZ())*0.5 << std::endl;
    break;	
  }
  out << "children [ " << std::endl;
  for(j=0;j<st->getPositions().size();j++)
    {
      apos = st->getPositions()[j];
      out << "Transform {  #" << apos->getVolume()->getName() << std::endl;
      switch(st->getStackType()){
      case sx:
	delta = apos->getBBX()/2 + deltap;
	out << "translation " << delta << " 0 " << " 0 " << std::endl;
	deltap += apos->getBBX();
	break;	
      case sy:
	delta = apos->getBBY()/2 + deltap;
	out << "translation " << " 0 " << delta << " 0 " << std::endl;
	deltap += apos->getBBY();
	break;	
      case sz:
	delta = apos->getBBZ()/2 + deltap;
	out << "translation " << " 0 " << " 0 " << delta << std::endl;
	deltap += apos->getBBZ();
	break;	
      }
      out << "children [ " << std::endl;
      apos->AcceptNotRec(this);
      out << "] " << std::endl; 
      out << "} " << std::endl;
    }  

  out << "] " << std::endl; 
  out << "} " << std::endl;
}



void  GDDVRMLSectionsVisitor::visitBox(GDDbox* box)
{
  out << "Shape {   #" << box->getName() << std::endl;

  string mat = box->getMaterial();

  if (mat == "CsI")
    {
      out << "appearance USE CSI" << std::endl;
    }
  else if (mat == "Si")
    {
      out << "appearance USE SI " << std::endl;
    }
  else if (mat == "rubber60A1")
    {
      out << "appearance USE RUBBER60A1 " << std::endl;
    }
  else if (mat == "FOAM05")
    {
      out << "appearance USE FOAM05 " << std::endl;
    }
  else if (mat == "Pb")
    {
      out << "appearance USE PB " << std::endl;
    }
  else if (mat == "C")
    {
      out << "appearance USE C " << std::endl;
    }
  else 
    {
      out << "appearance USE DUMMY " << std::endl;
    }
  out << "  geometry Box { " << std::endl;
  out << "                     size " 
	    << box->getX() << " " << box->getY() << " " << box->getZ() << std::endl; 
  out << "                    }" << std::endl;  
  out << "   }" << std::endl;
}

void  GDDVRMLSectionsVisitor::visitPosXYZ(GDDposXYZ* pos)
{
  out << "Transform { " << std::endl;
  out << "translation " << pos->getX() << " " <<
    pos->getY() << " " << pos->getZ() << 
    std::endl;
  out << "children [ " << std::endl;
  pos->getVolume()->AcceptNotRec(this);
  out << "] " << std::endl; 
  out << "} " << std::endl;
}


void  GDDVRMLSectionsVisitor::visitAxisPos(GDDaxisPos* pos)
{
  out << "Transform { " << std::endl;
  switch(pos->getPosDir()){
  case sx:
    out << "rotation " << " 1 0 0 " <<  pos->getRotation()*3.141927/180 << std::endl;  
    break;
  case sy:
    out << "rotation " << " 0 1 0 " <<  pos->getRotation()*3.141927/180 << std::endl;  
    break;
  case sz:
    out << "rotation " << " 0 0 1 " <<  pos->getRotation()*3.141927/180  << std::endl;  
    break;
  }
  out << "children [ " << std::endl;
  pos->getVolume()->AcceptNotRec(this);
  out << "] " << std::endl; 
  out << "} " << std::endl;

}

void  GDDVRMLSectionsVisitor::visitAxisMPos(GDDaxisMPos* pos)
{
  int i;
  int n;
  float g, w;


  n = pos->getNcopy();
  g = pos->getGap();
	  
  for(i=0;i<n;i++){
    out << "Transform { " << std::endl;
    switch(pos->getPosDir()){
    case sx:
      w = pos->getVolume()->getBBX();	
      out << "translation " << i*(w+g)-(n-1)*(w+g)*0.5  << " 0 " << " 0 " << std::endl;  
      break;
    case sy:
      w = pos->getVolume()->getBBY();
      out << "translation " << " 0 " << i*(w+g)-(n-1)*(w+g)*0.5  << " 0 " << std::endl;  
      break;
    case sz:
      w = pos->getVolume()->getBBZ();
      out << "translation " << " 0 " << " 0 " << i*(w+g)-(n-1)*(w+g)*0.5  << std::endl;  
      break;
    }
    out << "children [ " << std::endl;
    pos->getVolume()->AcceptNotRec(this);
    out << "] " << std::endl; 
    out << "} " << std::endl;
  }


}

void  GDDVRMLSectionsVisitor::visitIdField(GDDidField*)
{

}

void  GDDVRMLSectionsVisitor::visitSeg(GDDseg*)
{

}







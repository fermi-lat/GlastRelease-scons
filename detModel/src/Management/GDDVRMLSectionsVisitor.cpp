#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Management/GDDVRMLSectionsVisitor.h"
#include "detModel/Management/GDDmanager.h"
#include "detModel/Sections/GDDsection.h"
#include "detModel/Sections/GDDvolume.h"
#include "detModel/GDD.h"
#include "detModel/Sections/GDDbox.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDensemble.h"
#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDaxisMPos.h"
#include "detModel/Sections/GDDidField.h"
#include "detModel/Sections/GDDposition.h"
#include "detModel/Utilities/GDDcolorCreator.h"
#include "detModel/Utilities/GDDcolor.h"

/// This is the main constructor
GDDVRMLSectionsVisitor::GDDVRMLSectionsVisitor(std::string nvol)
{
  // Some variables
  unsigned int i;
  typedef std::map<std::string,GDDcolor*>M1;
  typedef std::map<std::string,int>M2;
  typedef std::map<std::string,GDDvolume*>M3;
  M3::const_iterator j;

  /// This visitor is not recursive
  setRecursive(0);
  actualVolume = nvol;
  depth = 0;

  /// This is the output file \todo The user should choose the name
  out.open("sections.wrl");
  
  /// We retrive the manager and the GDD
  GDDmanager* manager = GDDmanager::getPointer();
  GDD* g = manager->getGDD();
  
  /// We initialize the colors map for the material
  std::vector <std::string> names = g->getMaterialNames();
  GDDcolorCreator* cColor = new GDDcolorCreator(names.size());
  for(i=0;i<names.size();i++)
    colorsMap.insert(M1::value_type(names[i],cColor->getColor(i)));  
  
  /// We initialize the depth map
  M3 m = g->getVolumesMap();
  for(j=m.begin();j!=m.end();j++)
    depthMap.insert(M2::value_type(j->first,20));        
};

GDDVRMLSectionsVisitor::~GDDVRMLSectionsVisitor()
{
  out.close();
}

void GDDVRMLSectionsVisitor::visitGDD(GDD* gdd)
{
  typedef std::vector<GDDsection*>sec;
  std::vector <GDDsection*>::iterator i;
  
  out << "#VRML V2.0 utf8 " << std::endl;
  out << "#Generated by detModel " << std::endl;
  makeColor();
  
  sec s = gdd->getSections();
  for(i=s.begin(); i!=s.end();i++)
    (*i)->AcceptNotRec(this);
}

void  GDDVRMLSectionsVisitor::visitSection(GDDsection* section)
{
  float dimX, dimY, dimZ;
  GDDvolume* vol=0;
  std::vector <GDDvolume*>::iterator v;
  
  if (actualVolume == "")
    {
      std::vector<GDDvolume*> volumes = section->getVolumes();
      
      for(v=volumes.begin(); v!=volumes.end(); v++){
	if(GDDensemble* ens = dynamic_cast<GDDensemble*>(*v))
	  {
	    if (ens == section->getTopVolume())
	      vol = ens;
	  }
      }
    }
  else
    {
      GDDmanager* manager = GDDmanager::getPointer();
      if (manager->getGDD()->getVolumeByName(actualVolume))
	vol = manager->getGDD()->getVolumeByName(actualVolume);
      else
	{
	  std::cout << "No such volume" << std::endl;
	  exit(0);
	}
    }
  
  /// We calulate the dimensions for the projected views
  dimX = 2*vol->getBBox()->getXDim()/0.005;
  dimY = 2*vol->getBBox()->getYDim()/0.005;
  dimZ = 2*vol->getBBox()->getZDim()/0.005;

  /** 
   * We setup three standard view on the three coordinates planes;
   * since in VRML97 there is no orthographic camera, we simulate it
   * with a far away camera and a small field of view. \todo To be 
   * corrected; some views are not properly performed
   */
  
  out << "Viewpoint {" << std::endl;
  out << "position 0 0 " << dimX << std::endl;
  out << "orientation 0 0 1 0" << std::endl;
  out << "fieldOfView 0.005" << std::endl;
  out << "description \"Plane XY\"" << std::endl;
  out << "}" << std::endl;

  out << "Viewpoint {" << std::endl;
  out << "position 0 " << -dimY << " 0" << std::endl;
  out << "orientation 1 0 0 1.5707963" << std::endl;
  out << "fieldOfView 0.005" << std::endl;
  out << "description \"Plane XZ\"" << std::endl;
  out << "}" << std::endl;

  out << "Viewpoint {" << std::endl;
  out << "position " << -dimZ << " 0 0" << std::endl; 
  out << "orientation 0 1 0 -1.5707963" << std::endl;
  out << "fieldOfView 0.005" << std::endl;
  out << "description \"Plane ZY\"" << std::endl;
  out << "}" << std::endl;

  /// Here we start the visit of the hierarchy.
  vol->AcceptNotRec(this);
}


void  GDDVRMLSectionsVisitor::visitEnsemble(GDDensemble* ensemble)
{
  std::vector <GDDposition*>::iterator i;
  typedef std::vector<GDDposition*> pos;
  
  typedef std::map<std::string, int> M;
  M::const_iterator j; 

  j = depthMap.find(ensemble->getName());
  
  /// If we have not reached the depth limit we draw the full ensamble
  if(j->second > depth) 
    {
      depth++;
      
      out << "# " << ensemble->getName() << std::endl;
      
      /// Here the positioned volumes are visited
      pos p = ensemble->getPositions();
      for(i=p.begin(); i!=p.end();i++)
	(*i)->AcceptNotRec(this);
      
      /// Here the envelope is visited if the ensamble is a composition
      if (GDDcomposition* comp = dynamic_cast<GDDcomposition*>(ensemble))
	comp->getEnvelope()->AcceptNotRec(this);
      
      /// The depth level is decreased
      depth--;
    }
  else /// Else we draw only the bounding box
    {
      out << "Shape {   #" << ensemble->getName() << std::endl;
      out << "  geometry Box { " << std::endl;
      out << "                     size " 
	  << ensemble->getBBox()->getXDim() << " "
	  << ensemble->getBBox()->getXDim() << " "
	  << ensemble->getBBox()->getXDim() << std::endl;
      out << "                    }" << std::endl;  
      out << "   }" << std::endl;
    }
}

void  GDDVRMLSectionsVisitor::visitBox(GDDbox* box)
{
  out << "Shape {   #" << box->getName() << std::endl;
  out << "appearance USE " << box->getMaterial() << std::endl;
  
  out << "  geometry Box { " << std::endl;
  out << "                     size " 
      << box->getX() << " " << box->getY() << " " << box->getZ() << std::endl; 
  out << "                    }" << std::endl;  
  out << "   }" << std::endl;
}

void  GDDVRMLSectionsVisitor::visitPosXYZ(GDDposXYZ* pos)
{
  out << "Transform { " << std::endl;
  
  out << "translation " << pos->getX() << " " << pos->getY() << " " << pos->getZ() << std::endl;
  out << "rotation " << " 1 0 0 " <<  pos->getXRot()*3.141927/180 << std::endl;  
  out << "rotation " << " 0 1 0 " <<  pos->getYRot()*3.141927/180 << std::endl;  
  out << "rotation " << " 0 0 1 " <<  pos->getZRot()*3.141927/180 << std::endl;  
  out << "children [ " << std::endl;
  pos->getVolume()->AcceptNotRec(this);
  out << "] " << std::endl; 
  out << "} " << std::endl;
}


void  GDDVRMLSectionsVisitor::visitAxisMPos(GDDaxisMPos* pos)
{
  int i;
  int n;

  n = pos->getNcopy();

  for(i=0;i<n;i++){
    out << "Transform { " << std::endl;
    switch(pos->getAxisDir()){
    case (GDDstack::xDir):
      out << "rotation " << " 1 0 0 " <<  pos->getRotation()*3.141927/180 << std::endl;  
      out << "translation " << 
	pos->getDx()+pos->getDisp(i) << " " << pos->getDy() << " " << pos->getDz() << std::endl;
      break;
    case (GDDstack::yDir):
      out << "rotation " << " 0 1 0 " <<  pos->getRotation()*3.141927/180 << std::endl;  
      out << "translation " << 
	pos->getDx() << " " << pos->getDy()+pos->getDisp(i) << " " << pos->getDz() << std::endl;
      break;
    case (GDDstack::zDir):
      out << "rotation " << " 0 0 1 " <<  pos->getRotation()*3.141927/180 << std::endl;  
      out << "translation " << 
	pos->getDx() << " " << pos->getDy() << " " << pos->getDz()+pos->getDisp(i) << std::endl;
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
  /// We don't need identifiers in the VRML visitor
}


void  GDDVRMLSectionsVisitor::visitSeg(GDDseg*)
{
}

void GDDVRMLSectionsVisitor::setDepth(std::string name, int d)
{
  typedef std::map<std::string, int> M;
  M::const_iterator j; 
  
  j = depthMap.find(name);
  if (j == depthMap.end()) return;
  else 
    {
      depthMap.erase(name);
      depthMap.insert(M::value_type(name,d));
    }  
}

void GDDVRMLSectionsVisitor::setOpacity(std::string name, float op)
{
  typedef std::map<std::string, GDDcolor*> M;
  M::iterator j; 
  
  j = colorsMap.find(name);
  if (j == colorsMap.end()) return;
  else 
    j->second->setTra(op);
}

void GDDVRMLSectionsVisitor::setAllOpacity(float op)
{
  typedef std::map<std::string, GDDcolor*> M;
  M::const_iterator j; 
  
  for(j=colorsMap.begin();j!=colorsMap.end();j++)
    j->second->setTra(op);
}

void GDDVRMLSectionsVisitor::makeColor()
{
  typedef std::map<std::string, GDDcolor*> M;
  M::const_iterator j; 
  
  for(j=colorsMap.begin();j!=colorsMap.end();j++)
    {
      out << " DEF " << j->first << std::endl;
      out << " Appearance { " <<  std::endl;
      out << " material Material { " <<  std::endl;

      out << "       diffuseColor  " 
	  << j->second->getRed() << " " 
	  << j->second->getGreen() << " " 
	  << j->second->getBlue() 
	  << std::endl;
      out << "     transparency    " << j->second->getTra() << std::endl;          

      out << "     }" <<  std::endl;
      out << " }" <<  std::endl;
    }
}



#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#ifndef WIN32  // THB: not in our distribution
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#endif
#include "G4VSensitiveDetector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4SectionsVisitor.h"


#include "detModel/Management/Manager.h"
#include "detModel/Sections/Section.h"
#include "detModel/Sections/Choice.h"
#include "detModel/Gdd.h"
#include "detModel/Utilities/Vector.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Materials/Material.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Sections/IdField.h"

/// \todo Distruzione degli oggetti


G4SectionsVisitor::G4SectionsVisitor()
{
    typedef std::map<std::string,float>M1;
    typedef std::map<std::string,detModel::Material*>M2;
  
  M2::const_iterator k;

  setRecursive(0);
  setActualVolume("");

  detModel::Manager* manager = detModel::Manager::getPointer();
  detModel::Gdd* g = manager->getGdd();

  std::map<std::string, detModel::Material*> mat = g->getMaterials()->getMaterials();

  g->getMaterials()->generateColor();
  
  /// We initialize the opacity map
  for(k=mat.begin();k!=mat.end();k++)
    opacityMap.insert(M1::value_type(k->first,0.0));  

  compX = 0;
  compY = 0;
  compZ = 0;
  
  makeColor();
  halfPrec=0;
};

G4SectionsVisitor::~G4SectionsVisitor()
{
}


void G4SectionsVisitor::visitGdd(detModel::Gdd* gdd)
{
  (gdd->getSections() )[0]->AcceptNotRec(this);  
}

void  G4SectionsVisitor::visitSection(detModel::Section* section)
{
  detModel::Volume* topVolume = (section->getTopVolume());
  
  G4Material* ptMaterial = G4Material::GetMaterial("Galactic");

  G4Box* world
    = new G4Box(topVolume->getName(),
		2*topVolume->getBBox()->getXDim()*mm,
		2*topVolume->getBBox()->getYDim()*mm,
		2*topVolume->getBBox()->getZDim()*mm);
  
  worldlog = new G4LogicalVolume(world,ptMaterial,"motherVolume",0,0,0);

  /// We set the mother invisible
#ifndef WIN32 //THB
  worldlog->SetVisAttributes (G4VisAttributes::Invisible);  
#endif
  worldphys
    = new G4PVPlacement(0,G4ThreeVector(),"motherVolume",
			worldlog,0,false,0);


  if(actualVolume != "")
    {
      detModel::Gdd* gdd;
      detModel::Volume* vol;
      
      detModel::Manager* manager = detModel::Manager::getPointer();
      
      gdd = manager->getGdd();
      
      vol = gdd->getVolumeByName(actualVolume);
      if(vol!= 0)
	{
	  setActualMother(worldlog);
	  vol->AcceptNotRec(this);
	  
	  new G4PVPlacement(0,G4ThreeVector(),
			    g4Logicals[0],
			    g4Logicals[0]->GetName(),
			    worldlog,false,0);
	  std::cout<< "Geometry constructed " << std::endl;
	}
      else actualVolume = "";
    }
  if(actualVolume == "")
    {
      std::vector<detModel::Volume*>::iterator v;
      std::vector<detModel::Volume*> volumes = section->getVolumes();
      
      for(v=volumes.begin(); v!=volumes.end(); v++){
	  if(detModel::Ensemble* ens = dynamic_cast<detModel::Ensemble*>(*v))
	    {
	      if (ens == section->getTopVolume())
		{
		  setActualMother(worldlog);
		  ens->AcceptNotRec(this);
	      
		  new G4PVPlacement(0,G4ThreeVector(),
				    g4Logicals[0],
				    g4Logicals[0]->GetName(),
				    worldlog,false,0);
		}
	    }
      }
      std::cout<< "Geometry constructed " << std::endl;
    }
  
  
}

void  G4SectionsVisitor::visitEnsemble(detModel::Ensemble* ensemble)
{
  unsigned int i, logind;
  G4LogicalVolume *temp=actualMother;

  G4Material* ptMaterial = G4Material::GetMaterial("Galactic");

  g4Boxes.push_back(
		    new G4Box(ensemble->getName(),
			      ensemble->getBBox()->getXDim()*mm/2,
			      ensemble->getBBox()->getYDim()*mm/2,
			      ensemble->getBBox()->getZDim()*mm/2)
		    );


  g4Logicals.push_back(
		       new G4LogicalVolume(
					   g4Boxes.back(),ptMaterial,ensemble->getName(),
					   0,0,0)
		       );


  /// We don't want to see the composition volume
#ifndef WIN32 //THB
  g4Logicals.back()->SetVisAttributes (G4VisAttributes::Invisible);
#endif
  /// This is the index of the actual physic volume
  logind = g4Logicals.size() - 1;

  for(i=0; i<ensemble->getPositions().size();i++){
    detModel::Position* pos = ensemble->getPositions()[i];
    setActualMother(g4Logicals[logind]);
    pos->AcceptNotRec(this);
  }


  actualMother=temp;
}


void  G4SectionsVisitor::visitBox(detModel::Box* box)
{
    typedef std::map<std::string, G4VisAttributes*> M;
  M::const_iterator i; 

  G4Material* ptMaterial = G4Material::GetMaterial(box->getMaterial());
  g4Boxes.push_back(
		    new G4Box(box->getName(),
			      box->getX()*mm/2,
			      box->getY()*mm/2,
			      box->getZ()*mm/2)
		    );

  g4Logicals.push_back(
		       new G4LogicalVolume(
					   g4Boxes.back(),ptMaterial,box->getName(),
					   0,0,0)
		       );



  i = g4VisAttributes.find(box->getMaterial());
  if(i != g4VisAttributes.end()) 
    g4Logicals.back()->SetVisAttributes(i->second);
  
}  


void  G4SectionsVisitor::visitPosXYZ(detModel::PosXYZ* pos)
{
  unsigned int i;
  int  ind = -1;
  char temp[10];  std::string tmp = "";
  typedef std::map<G4VPhysicalVolume*,std::string>M;
  std::string volName;

  if (detModel::Choice* choice = dynamic_cast<detModel::Choice*>(pos->getVolume()))
    {
      detModel::Manager* man = detModel::Manager::getPointer();
      volName = (choice->getVolumeByMode(man->getMode()))->getName();
    }
  else
    volName = pos->getVolume()->getName();

  
  for(i=0;i<g4Logicals.size();i++){
    if (g4Logicals[i]->GetName() == volName) 
      {
	ind = i;
	break;
      }
  }//endfor
  if (ind == -1)
    {
      ind = g4Logicals.size();
      pos->getVolume()->AcceptNotRec(this);
    }

  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(((pos->getXRot())*GDDPI)/180.0);
  rm->rotateY(((pos->getYRot())*GDDPI)/180.0);
  rm->rotateZ(((pos->getZRot())*GDDPI)/180.0);

  g4Physicals.push_back( new  G4PVPlacement(rm, G4ThreeVector((pos->getX())*mm, 
							      (pos->getY())*mm, 
							      (pos->getZ())*mm),
					    g4Logicals[ind],
					    volName,
					    actualMother,false,0));
  
  for(i=0;i<(pos->getIdFields()).size();i++)
    { 
      sprintf(temp, "/%d", (int) pos->getIdFields()[i]->getValue()); 
      tmp = tmp+temp;
    }

  if(pos->getIdFields().size() != 0)
    g4Identifiers.insert(M::value_type(g4Physicals.back(),tmp));  

}


void  G4SectionsVisitor::visitAxisMPos(detModel::AxisMPos* pos)
{
  unsigned int i, ncopy, j;
  int  ind = -1;
  std::string tmp = "";
  char temp[10];
  typedef std::map<G4VPhysicalVolume*,std::string>M;
  detModel::Volume* volume;


  if (detModel::Choice* choice = dynamic_cast<detModel::Choice*>(pos->getVolume()))
    {
      detModel::Manager* man = detModel::Manager::getPointer();
      volume = choice->getVolumeByMode(man->getMode());
    }
  else
    volume = pos->getVolume();

  std::string volName = volume->getName();
 

  ncopy = pos->getNcopy();

  for(i=0;i<g4Logicals.size();i++)
    if (g4Logicals[i]->GetName() == volName) 
      {
	ind = i;
	break;
      }
  
  if (ind == -1)
    {
      ind = g4Logicals.size();
      pos->getVolume()->AcceptNotRec(this);
    }

  G4RotationMatrix* rm = new G4RotationMatrix();

  for(i=0;i<ncopy;i++)
    {
  
      switch(pos->getAxisDir()){
      case (detModel::Stack::xDir):
	{
	  rm->rotateX((pos->getRotation()*GDDPI)/180.0);
	  g4Physicals.push_back(
				new G4PVPlacement(rm, G4ThreeVector((pos->getDx()+pos->getDisp(i))*mm, 
								    (pos->getDy())*mm, 
								    (pos->getDz())*mm),
						  getLogicalByName(volName),
						  volName,
						  actualMother,false,0)
				);	
	}
	break;
      case (detModel::Stack::yDir):
	{
	  rm->rotateY(pos->getRotation()*GDDPI/180.0);
	  g4Physicals.push_back(
				new G4PVPlacement(rm, G4ThreeVector((pos->getDx())*mm, 
								    (pos->getDy()+pos->getDisp(i))*mm, 
								    (pos->getDz())*mm),
						  getLogicalByName(volName),
						  volName,
						  actualMother,false,0)
				);
	}	
	break;
      case (detModel::Stack::zDir):
	{
	  rm->rotateZ(pos->getRotation()*GDDPI/180.0);
	  g4Physicals.push_back(
				new G4PVPlacement(rm, G4ThreeVector((pos->getDx())*mm, 
								    (pos->getDy())*mm, 
								    (pos->getDz()+pos->getDisp(i))*mm),
						  getLogicalByName(volName),
						  volName,
						  actualMother,false,0)
				);
	}	
	break;
      }//endSwitch
    }//endfor


  for(j=0;j<pos->getIdFields().size();j++)
    {	
      sprintf(temp, "/%d", 
	      (int)(pos->getIdFields()[j]->getValue())+(int)(pos->getIdFields()[j]->getStep()*i)); 
      tmp = tmp+temp;
    }
  if(pos->getIdFields().size() != 0)
    g4Identifiers.insert(M::value_type(g4Physicals.back(),tmp));  
}

void  G4SectionsVisitor::visitIdField(detModel::IdField*)
{

}

void  G4SectionsVisitor::visitSeg(detModel::Seg*)
{

}

void G4SectionsVisitor::setOpacity(std::string name, float op)
{
    typedef std::map<std::string, float> M;
  M::const_iterator j; 
  
  j = opacityMap.find(name);
  if (j == opacityMap.end()) return;
  else 
    {
      opacityMap.erase(name);
      opacityMap.insert(M::value_type(name,op));
    }  
}

void G4SectionsVisitor::makeColor()
{
#ifndef WIN32  //THB
    typedef std::map<std::string, G4VisAttributes*> Ma;
  Ma::const_iterator j; 
  typedef std::map<std::string, detModel::Color*> M;
  M::const_iterator k; 
  
  detModel::Manager* manager = detModel::Manager::getPointer();
  detModel::Gdd* gdd = manager->getGdd();

  
  std::map<std::string,detModel::Color*> colorsMap = gdd->getMaterials()->getMaterialColors();
  
  for(k=colorsMap.begin();k!=colorsMap.end();k++)
    {
      g4VisAttributes.insert(Ma::value_type(k->first,
					    new G4VisAttributes(G4Colour(
									 (k->second)->getRed(),
									 (k->second)->getGreen(),
									 (k->second)->getBlue()))));
    }

  for (j=g4VisAttributes.begin(); j!=g4VisAttributes.end(); j++){
    j->second->SetForceSolid(true);
  }
#endif
}

G4LogicalVolume* G4SectionsVisitor::getLogicalByName(std::string name){
  unsigned int i;

  for(i=0;i<g4Logicals.size();i++)
    if (g4Logicals[i]->GetName() == name) return g4Logicals[i];

  return 0;
}










#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <iomanip>

#include "detModel/Management/Manager.h"
#include "detModel/Sections/Section.h"
#include "detModel/Gdd.h"
#include "detModel/Sections/Choice.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Tube.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Sections/IdField.h"
#include "detModel/Sections/Position.h"
#include "detModel/Utilities/PositionedVolume.h"
#include "detModel/Utilities/Vector.h"

#include "detModel/Management/IDmapBuilder.h"
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"


namespace detModel{

  IDmapBuilder::IDmapBuilder(std::string nvol)
  {
    setRecursive(0);
    m_actualVolume = nvol;
  };

IDmapBuilder::~IDmapBuilder()
{
  for(std::map<idents::VolumeIdentifier,const PositionedVolume*>::iterator i=m_volMap.begin();
      i != m_volMap.end();i++)
    {
      delete (*i).second;
      (*i).second = 0;
    }
	
  m_volMap.clear();
}


void IDmapBuilder::visitGdd(Gdd* gdd)
{
  typedef std::vector<Section*>sec;
  std::vector <Section*>::iterator i;

  m_gdd = gdd;
  sec s = gdd->getSections();
  for(i=s.begin(); i!=s.end();i++)
    (*i)->AcceptNotRec(this);

}
  
void  IDmapBuilder::visitSection(Section* section)
{
  Volume* vol=0;
  std::vector <Volume*>::iterator v;
  
  if (m_actualVolume == "")
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
      //      Manager* manager = Manager::getPointer();
      //      if (manager->getGdd()->getVolumeByName(m_actualVolume))
      //	vol = manager->getGdd()->getVolumeByName(m_actualVolume);
      if (m_gdd->getVolumeByName(m_actualVolume))
	vol = m_gdd->getVolumeByName(m_actualVolume);
      else
	{
	  std::cout << "No such volume" << std::endl;
	  exit(0);
	}
    }

  vol->AcceptNotRec(this);
}

void  IDmapBuilder::visitEnsemble(Ensemble* ensemble)
{
  std::vector <Position*>::iterator i;
  typedef std::vector<Position*> pos;
  
  /// Here the positioned volumes are visited
  pos p = ensemble->getPositions();
  for(i=p.begin(); i!=p.end();i++)
    (*i)->AcceptNotRec(this);

}


void  IDmapBuilder::visitBox(Box*)
{
}

void  IDmapBuilder::visitTube(Tube*)
{
}

void  IDmapBuilder::visitPosXYZ(PosXYZ* pos)
{
  HepRotation tempRotation = m_actualRot;
  Hep3Vector tempPos = m_actualPos;
  unsigned int i;
  idents::VolumeIdentifier tempID = m_actualID;

  m_actualPos = m_actualPos + m_actualRot*Hep3Vector(pos->getX(), pos->getY(), pos->getZ());

  m_actualRot.rotateX(pos->getXRot()*GDDPI/180);
  m_actualRot.rotateY(pos->getYRot()*GDDPI/180);
  m_actualRot.rotateZ(pos->getZRot()*GDDPI/180);
  
  /// Set the identifier
  for(i=0;i<(pos->getIdFields()).size();i++)
    m_actualID.append((int) pos->getIdFields()[i]->getValue());
  
  insertVolume(pos->getVolume());
  pos->getVolume()->AcceptNotRec(this);
  
  m_actualID = tempID;

  m_actualPos = tempPos;
  m_actualRot = tempRotation;
}


void  IDmapBuilder::visitAxisMPos(AxisMPos* pos)
{
  unsigned int i, j;
  HepRotation tempRotation = m_actualRot;
  Hep3Vector tempPos = m_actualPos;
  idents::VolumeIdentifier tempID = m_actualID;

  Hep3Vector stackDir; // unit vector to set 
  
  switch(pos->getAxisDir()){
  case (detModel::Stack::xDir):
    {    
      m_actualRot.rotateX(-pos->getRotation()*GDDPI/180); 
      stackDir=Hep3Vector(1,0,0); 
    }
      break;
  case (detModel::Stack::yDir):
    {
      m_actualRot.rotateY(-pos->getRotation()*GDDPI/180); 
      stackDir=Hep3Vector(0,1,0); 
    }
    break;
  case (detModel::Stack::zDir):
    {
      m_actualRot.rotateZ(-pos->getRotation()*GDDPI/180); 
      stackDir=Hep3Vector(0,0,1); 
    }
    break;
  }

  Hep3Vector origin(pos->getDx(),pos->getDy(),pos->getDz());
  

  for(i=0;i<pos->getNcopy();i++)
    {
      m_actualPos = m_actualPos + origin + (pos->getDisp(i))*(m_actualRot*stackDir);

      /// ID stuff
      for(j=0;j<pos->getIdFields().size();j++)
	{	
	  m_actualID.append((int)(pos->getIdFields()[j]->getValue())+(int)(pos->getIdFields()[j]->getStep()*i));
	}


      insertVolume(pos->getVolume());      
      pos->getVolume()->AcceptNotRec(this);
      
      m_actualID = tempID;
      m_actualPos = tempPos;
    }
  m_actualRot= tempRotation;
}

void  IDmapBuilder::visitIdField(IdField*)
{

}

void  IDmapBuilder::visitSeg(Seg*)
{

}

void IDmapBuilder::insertVolume(Volume* vol)
{
  Volume* volume;
  /// This is necessary to deal with Choice volume
  if(Choice* choice = dynamic_cast<Choice*>(vol))
    {
      //      Manager* man = Manager::getPointer();
      //      volume = choice->getVolumeByMode(man->getMode());
      volume = choice->getVolumeByMode(m_mode);
    }
  else volume = vol;

  /// This is the criterium to add the volume
  if(Shape* shape = dynamic_cast<Shape*>(volume))
    {
      if (shape->getSensitive())
	{    
	  PositionedVolume* pos = new PositionedVolume();
	  pos->setVolume(volume);
	  pos->setTranslation(m_actualPos);
	  pos->setRotation(m_actualRot);
	  m_volMap[m_actualID] = pos;
          // also save a vector of id's
          m_idvec.push_back(m_actualID);
	}
    }
}  

const std::map<idents::VolumeIdentifier, const PositionedVolume*>* IDmapBuilder::getVolMap() const
{
  return &m_volMap;
}

const PositionedVolume* IDmapBuilder::getPositionedVolumeByID(idents::VolumeIdentifier id) const
{
    PVmap::const_iterator i = m_volMap.find(id);
    return i != end() ?  (*i).second : 0;

}

void IDmapBuilder::summary(std::ostream & out) 
{
    std::map<std::string, unsigned int> names;

    for( PVmap::const_iterator it = begin(); it !=end(); ++it){
        const detModel::PositionedVolume * pv = (*it).second;
        const detModel::Volume* vol = pv->getVolume();
        std::string name = vol->getName();
        names[name]++;
    }
    out << "Summary of IDmap contents:" << std::endl
        << std::setw(20) << "Detector name" << std::setw(10) << "count" << std::endl
        << std::setw(20) << "-------------" << std::setw(10) << "-----" << std::endl;
    
    for( std::map<std::string, unsigned int>::iterator j=names.begin(); j!=names.end(); ++j){
        out << std::setw(20) << (*j).first << std::setw(10) << (*j).second << std::endl;
    }
    out << std::setw(20) << "total:" << std::setw(10) << size() << std::endl;
}

IDmapBuilder::PVmap::const_iterator IDmapBuilder::begin()const {
  return m_volMap.begin();
}
IDmapBuilder::PVmap::const_iterator IDmapBuilder::end()const {
  return m_volMap.end();
}

size_t IDmapBuilder::size()const {
  return m_volMap.size();
}

}//namespace

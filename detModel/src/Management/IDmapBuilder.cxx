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
#include "detModel/Sections/Volume.h"
#include "detModel/Sections/Shape.h"

#include "detModel/Utilities/PositionedVolume.h"
#include "detModel/Utilities/Vector.h"

#include "detModel/Management/IDmapBuilder.h"
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"


namespace detModel{

  IDmapBuilder::IDmapBuilder(std::string nvol) 
    : m_topSeen(false), 
    m_transformPrefix(HepTransform3D::Identity),
    m_inverseTransformPrefix(HepTransform3D::Identity)
  {
    setRecursive(0);
    m_actualVolume = nvol;
  };

IDmapBuilder::~IDmapBuilder()
{
  for (std::map<idents::VolumeIdentifier,const PositionedVolume*>::iterator i
         = m_volMap.begin(); i != m_volMap.end(); i++)
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
  for (i=s.begin(); i!=s.end();i++)
    (*i)->AcceptNotRec(this);

}
  
void  IDmapBuilder::visitSection(Section* section)
{
  Volume* vol=0;
  std::vector <Volume*>::iterator v;

  if (m_actualVolume != "") {  // just make sure it exists
    if (!(m_gdd->getVolumeByName(m_actualVolume))) {
      std::cout << "No such volume" << std::endl;
      exit(0);
    }
  }
  //  if (m_actualVolume == "")
  // Now always start with section's top volume
  std::vector<Volume*> volumes = section->getVolumes();
      
  for (v=volumes.begin(); v!=volumes.end(); v++){
    if (Ensemble* ens = dynamic_cast<Ensemble*>(*v))
    {
      if (ens == section->getTopVolume()) {
        vol = ens;
        break;         // I think! must ask Riccardo if break is OK here
      }
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
  for (i=p.begin(); i!=p.end();i++)
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

  m_actualPos = m_actualPos + 
    m_actualRot*Hep3Vector(pos->getX(), pos->getY(), pos->getZ());

  m_actualRot.rotateX(pos->getXRot()*GDDPI/180);
  m_actualRot.rotateY(pos->getYRot()*GDDPI/180);
  m_actualRot.rotateZ(pos->getZRot()*GDDPI/180);
  
  /// Set the identifier
  for (i=0;i<(pos->getIdFields()).size();i++)
    m_actualID.append((int) pos->getIdFields()[i]->getValue());
  
  Volume* vol = pos->getVolume();
  if (!(m_topSeen)) {
    if (vol->getName() == m_actualVolume) {
      m_topSeen = true;
      m_IDPrefix = tempID;
      m_transformPrefix = HepTransform3D(tempRotation, tempPos);
      m_inverseTransformPrefix = m_transformPrefix.inverse();
    }
  }

  insertVolume(vol);
  vol->AcceptNotRec(this);
  
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
      m_actualRot.rotateX(pos->getRotation()*GDDPI/180); 
      stackDir=Hep3Vector(1,0,0); 
    }
      break;
  case (detModel::Stack::yDir):
    {
      m_actualRot.rotateY(pos->getRotation()*GDDPI/180); 
      stackDir=Hep3Vector(0,1,0); 
    }
    break;
  case (detModel::Stack::zDir):
    {
      m_actualRot.rotateZ(pos->getRotation()*GDDPI/180); 
      stackDir=Hep3Vector(0,0,1); 
    }
    break;
  }

  Hep3Vector origin(pos->getDx(),pos->getDy(),pos->getDz());
  

  for(i=0;i<pos->getNcopy();i++)
    {
      m_actualPos = m_actualPos + origin + 
        (pos->getDisp(i))*(m_actualRot*stackDir);

      /// ID stuff
      for(j=0;j<pos->getIdFields().size();j++)
	{	
	  m_actualID.append((int)(pos->getIdFields()[j]->getValue())+
                            (int)(pos->getIdFields()[j]->getStep()*i));
	}

      Volume* vol = pos->getVolume();
      if (!(m_topSeen)) {
        if (vol->getName() == m_actualVolume) {
          m_topSeen = true;
          m_IDPrefix = tempID;
          m_transformPrefix = HepTransform3D(tempRotation, tempPos);
          m_inverseTransformPrefix = m_transformPrefix.inverse();
        }
      }

      insertVolume(vol);
      vol->AcceptNotRec(this);
      
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

  /// This is the criterion to add the volume
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

const std::map<idents::VolumeIdentifier, const PositionedVolume*>* 
IDmapBuilder::getVolMap() const
{
  return &m_volMap;
}

const PositionedVolume* 
IDmapBuilder::getPositionedVolumeByID(idents::VolumeIdentifier id) const
{
    PVmap::const_iterator i = m_volMap.find(id);
    return i != end() ?  (*i).second : 0;

}



bool IDmapBuilder::getTransform3DByID(idents::VolumeIdentifier id, 
                                      HepTransform3D* tr)
{
  const PositionedVolume* pv = getPositionedVolumeByID(id);
  if (pv)
    {
      HepTransform3D temp(pv->getRotation(), pv->getTranslation());
      *tr = temp;
      return true;
    }
  else return false;
}

bool IDmapBuilder::getTransform3DByTopID(idents::VolumeIdentifier id, 
                                         HepTransform3D* tr)
{
  idents::VolumeIdentifier fullId = m_IDPrefix;
  fullId.append(id);
  return getTransform3DByID(fullId, tr);
}

bool IDmapBuilder::getTopTransform3DByTopID(idents::VolumeIdentifier id, 
                                            HepTransform3D* tr) {

  // Don't do unnecessary algebra for common case of top volume = world
  if (m_IDPrefix.size() == 0) return getTransform3DByID(id, tr);

  if (!(getTransform3DByTopID(id, tr))) return false;
  *tr = m_inverseTransformPrefix * (*tr);

  return true;
}

bool IDmapBuilder::getTopTransform3DByID(idents::VolumeIdentifier id, 
                                         HepTransform3D* tr) {
  // Don't do unnecessary algebra for common case of top volume = world
  if (m_IDPrefix.size() == 0) return getTransform3DByID(id, tr);

  if (!(getTransform3DByID(id, tr))) return false;
  *tr = m_inverseTransformPrefix * (*tr);

  return true;
}
  

bool IDmapBuilder::getShapeByID(idents::VolumeIdentifier id,
                                std::string* s, std::vector<double>* params)
{
  const PositionedVolume* pv = getPositionedVolumeByID(id);
  if (pv)
    {
      Volume* volume = pv->getVolume();
      if (Box* box = dynamic_cast<Box*>(volume))
        {
          *s = "box";
          params->push_back(box->getX());
          params->push_back(box->getY());
          params->push_back(box->getZ());
        }
      else
        return false;
      return true;
    }
  else return false;
}

bool IDmapBuilder::getShapeByTopID(idents::VolumeIdentifier relId, 
                                   std::string* s, 
                                   std::vector<double>* params) {

  if (m_IDPrefix.size() == 0) return getShapeByID(relId, s, params);

  idents::VolumeIdentifier global = m_IDPrefix;
  global.append(relId);

  return getShapeByID(global, s, params);
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
        << std::setw(20) << "Detector name";
    out << std::setw(10) << "count" << std::endl
        << std::setw(20) << "-------------" 
        << std::setw(10) << "-----" << std::endl;
    
    for( std::map<std::string, unsigned int>::iterator j=names.begin(); 
         j!=names.end(); ++j){
      out << std::setw(20) << (*j).first;
      out << std::setw(10) << (*j).second << std::endl;
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

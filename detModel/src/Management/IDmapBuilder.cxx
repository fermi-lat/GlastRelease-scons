#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

#include "detModel/Management/Manager.h"
#include "detModel/Sections/Section.h"
#include "detModel/Gdd.h"
#include "detModel/Sections/Choice.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Tube.h"
#include "detModel/Sections/Trap.h"
#include "detModel/Sections/Sphere.h"
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
  using idents::VolumeIdentifier;

  // used below in orderSensitive
  typedef  std::pair<IDmapBuilder::PVmap::const_iterator, float> MapIt_Float;
  
  // Function object used to sort MapIt_Float
  class MapItCompare {
  public:
    bool operator() (const MapIt_Float& a, const MapIt_Float& b) {
      return a.second < b.second;
    }
  };

  IDmapBuilder::IDmapBuilder(std::string nvol) 
    : m_topSeen(false), 
    m_transformPrefix(HepGeom::Transform3D::Identity),
    m_inverseTransformPrefix(HepGeom::Transform3D::Identity)
  {
    setRecursive(0);
    m_actualVolume = nvol;
  };

  IDmapBuilder::~IDmapBuilder()
  {
    for (std::map<VolumeIdentifier,const PositionedVolume*>::iterator i
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

  void  IDmapBuilder::visitTrap(Trap*)
  {
  }

  void  IDmapBuilder::visitSphere(Sphere*)
  {
  }

  void  IDmapBuilder::visitPosXYZ(PosXYZ* pos)
  {
    CLHEP::HepRotation tempRotation = m_actualRot;
    CLHEP::Hep3Vector tempPos = m_actualPos;
    unsigned int i;
    VolumeIdentifier tempID = m_actualID;

    m_actualPos = m_actualPos + 
      m_actualRot*CLHEP::Hep3Vector(pos->getX(), pos->getY(), pos->getZ());

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
        //      m_IDPrefix = tempID;
        //      m_transformPrefix = HepTransform3D(tempRotation, tempPos);
        //  Try including identifier(s) & transform for top volume as
        //  part of "prefix"
        m_IDPrefix = m_actualID;
        m_transformPrefix = HepGeom::Transform3D(m_actualRot, m_actualPos);
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
    CLHEP::HepRotation tempRotation = m_actualRot;
    CLHEP::Hep3Vector tempPos = m_actualPos;
    VolumeIdentifier tempID = m_actualID;

    CLHEP::Hep3Vector stackDir; // unit vector to set 
  
    switch(pos->getAxisDir()){
    case (detModel::Stack::xDir):
      {    
        m_actualRot.rotateX(pos->getRotation()*GDDPI/180); 
        stackDir=CLHEP::Hep3Vector(1,0,0); 
      }
      break;
    case (detModel::Stack::yDir):
      {
        m_actualRot.rotateY(pos->getRotation()*GDDPI/180); 
        stackDir=CLHEP::Hep3Vector(0,1,0); 
      }
      break;
    case (detModel::Stack::zDir):
      {
        m_actualRot.rotateZ(pos->getRotation()*GDDPI/180); 
        stackDir=CLHEP::Hep3Vector(0,0,1); 
      }
      break;
    }

    CLHEP::Hep3Vector origin(pos->getDx(),pos->getDy(),pos->getDz());
  

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
          m_IDPrefix = m_actualID;
          m_transformPrefix = HepGeom::Transform3D(m_actualRot, m_actualPos);
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
    else if (Composition* cmp = dynamic_cast<Composition*>(volume)) {
      Shape* shape = cmp->getEnvelope();

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

  const std::map<VolumeIdentifier, const PositionedVolume*>* 
  IDmapBuilder::getVolMap() const
  {
    return &m_volMap;
  }

  const PositionedVolume* 
  IDmapBuilder::getPositionedVolumeByID(VolumeIdentifier id) const
  {
    PVmap::const_iterator i = m_volMap.find(id);
    return i != end() ?  (*i).second : 0;

  }



  bool IDmapBuilder::getTransform3DByID(VolumeIdentifier id, 
                                        HepGeom::Transform3D* tr)
  {
    const PositionedVolume* pv = getPositionedVolumeByID(id);
    if (pv)
    {
      HepGeom::Transform3D temp(pv->getRotation(), pv->getTranslation());
      *tr = temp;
      return true;
    }
    else return false;
  }

  bool IDmapBuilder::getTransform3DByTopID(VolumeIdentifier id, 
                                           HepGeom::Transform3D* tr)
  {
    VolumeIdentifier fullId = m_IDPrefix;
    fullId.append(id);
    return getTransform3DByID(fullId, tr);
  }

  bool IDmapBuilder::getTopTransform3DByTopID(VolumeIdentifier id, 
                                              HepGeom::Transform3D* tr) {

    // Don't do unnecessary algebra for common case of top volume = world
    if (m_IDPrefix.size() == 0) return getTransform3DByID(id, tr);

    if (!(getTransform3DByTopID(id, tr))) return false;
    *tr = m_inverseTransformPrefix * (*tr);

    return true;
  }

  bool IDmapBuilder::getTopTransform3DByID(VolumeIdentifier id, 
                                           HepGeom::Transform3D* tr) {
    // Don't do unnecessary algebra for common case of top volume = world
    if (m_IDPrefix.size() == 0) return getTransform3DByID(id, tr);

    if (!(getTransform3DByID(id, tr))) return false;
    *tr = m_inverseTransformPrefix * (*tr);

    return true;
  }
  

  bool IDmapBuilder::getShapeByID(VolumeIdentifier id,
                                  std::string* s, std::vector<double>* params)
  {
    const PositionedVolume* pv = getPositionedVolumeByID(id);
    if (pv)
    {
      Volume* volume = pv->getVolume();
      if  (Composition* cmp = dynamic_cast<Composition*>(volume)) {
        volume = cmp->getEnvelope();
      }
        
      if (Box* box = dynamic_cast<Box*>(volume))
      {
        *s = "box";
        params->push_back(box->getX());
        params->push_back(box->getY());
        params->push_back(box->getZ());
      }
      else if (Tube* tube = dynamic_cast<Tube*>(volume) )
      {
        *s="tube";
        params->push_back(tube->getZ());
        params->push_back(tube->getRin());
        params->push_back(tube->getRout());
      }
      else if (Trap* trap = dynamic_cast<Trap*>(volume) )
      {
        *s="trap";
        params->push_back(trap->getX1());    // get all X-related params first
        params->push_back(trap->getX2());
        params->push_back(trap->getXDiff());
        params->push_back(trap->getY());
        params->push_back(trap->getZ());
      }
      else if (Sphere* sphere = dynamic_cast<Sphere*>(volume) )
      {
        *s="sphere";
        params->push_back(sphere->getRin());
        params->push_back(sphere->getRout());
        params->push_back(sphere->getPhiMin());
        params->push_back(sphere->getPhiMax());
        params->push_back(sphere->getThetaMin());
        params->push_back(sphere->getThetaMax());
      }      
      return true;
    }
    else return false;
  }

  bool IDmapBuilder::getShapeByTopID(VolumeIdentifier relId, 
                                     std::string* s, 
                                     std::vector<double>* params) {

    if (m_IDPrefix.size() == 0) return getShapeByID(relId, s, params);

    VolumeIdentifier global = m_IDPrefix;
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

  unsigned IDmapBuilder::orderSensitive(std::vector<VolumeIdentifier>& vols,
                                        const VolumeIdentifier& sample,
                                        const std::vector<bool> wild,
                                        axisDir coord, 
                                        bool increasing) const {
    unsigned idSize = (unsigned) sample.size();
    if (idSize != wild.size() ) { // error
      std::string 
        errtxt("IDmapBuilder::orderSensitive: incompatible arg sizes");
      throw std::invalid_argument(errtxt);
    }
    // Build mask.  Will be used to find set of VolumeIdentifiers
    // satisfying the constraints.  
    // Note: this only works because all identifiers are packed into 
    // an int64 the same way.
    VolumeIdentifier mask;
    unsigned maxFieldVal = VolumeIdentifier::maxFieldValue();
    for (unsigned ix = 0; ix < idSize; ix ++) {
      if (wild[ix]) mask.append(0);
      else mask.append(maxFieldVal);
    }
    VolumeIdentifier::int64 maskValue = mask.getValue();
    VolumeIdentifier maskedSample;  
    maskedSample.init(maskValue & sample.getValue(), sample.size());
    

    // Scan through map.  Save iterators and assoc. coord of interest.

    std::vector<MapIt_Float> matches;
    matches.reserve(10);

    PVmap::const_iterator mapIt = begin();

    while (mapIt != end() ) {
      if (mapIt->first.size() == sample.size() ) {
        if ((mapIt->first.getValue() & maskValue) == 
            maskedSample.getValue() ) {   // got one
          CLHEP::Hep3Vector trans = mapIt->second->getTranslation();
          float coordVal;
          switch(coord) {
          case xDir: 
            coordVal = trans.x();
            break;
          case yDir: 
            coordVal = trans.y();
            break;
          case zDir: 
            coordVal = trans.z();
            break;
          }

          MapIt_Float match(mapIt, coordVal);
          matches.push_back(match);
        }
      }
      ++mapIt;
    }
    MapItCompare cmp;
    std::sort(matches.begin(), matches.end(), cmp);

    if (increasing) {
      for (unsigned ix = 0; ix < matches.size(); ix++) {
        vols.push_back((matches[ix].first)->first);
      }
    }
    else {
      for (unsigned ix = matches.size(); ix > 0; ix--) { 
        vols.push_back(matches[ix-1].first->first);
      }
    }
      
    return 1;
  }

}//namespace

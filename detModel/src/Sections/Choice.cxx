#include <string>
#include <map>
#include "detModel/Sections/Choice.h"
#include "detModel/Sections/BoundingBox.h"
#include "detModel/Sections/Boundable.h"
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Management/Manager.h"

namespace detModel{

  void Choice::addCaseName(std::string name, std::string volume)
  {
    typedef std::map<std::string, std::string> M;
    
    casesName.insert(M::value_type(name, volume));
  }

  std::string Choice::getVolumeNameByMode(std::string mode)
  {
    typedef std::map<std::string, std::string> M;
    M::const_iterator i = casesName.find(mode);
    if(i == casesName.end()) return defaultName;
    else return i->second;
  }

  void Choice::addCase(std::string name, Volume* volume)
  {
    typedef std::map<std::string, Volume* > M;
    
    cases.insert(M::value_type(name, volume));
  }

  Volume* Choice::getVolumeByMode(std::string mode)
  {
    typedef std::map<std::string, Volume*> M;
    M::const_iterator i; 
    
    i = cases.find(mode);
    if(i == cases.end()) return defaultVol;
    else return i->second;
  }

  void Choice::Accept(SectionsVisitor* s)
  {
    //    Manager* man = Manager::getPointer();
    //    Volume* v = getVolumeByMode(man->getMode());
    Volume* v = getVolumeByMode(s->getMode());
    v->Accept(s);
  }
  
  void Choice::AcceptNotRec(SectionsVisitor* s)
  {
    //    Manager* man = Manager::getPointer();
    //    Volume* v = getVolumeByMode(man->getMode());
    Volume* v = getVolumeByMode(s->getMode());
    v->AcceptNotRec(s);

  }  

  void Choice::buildBB()
  {
    BoundingBox* bb = getBBox();
    bb->merge(defaultVol->getBBox(),0,0,0);
  }
  
}





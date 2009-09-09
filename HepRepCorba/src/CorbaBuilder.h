#ifndef CORBABUILDER_H
#define CORBABUILDER_H
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include "HepRepSvc/IBuilder.h"

#include "HepRep.hh"

/** 
 *  @class CorbaBuilder
 *
 *  @brief 
 *
 *  @author 
 */

class CorbaBuilder : public IBuilder{
 public:
  //  CorbaBuilder();

  virtual ~CorbaBuilder();

  virtual void addAttDef(std::string name,
                           std::string desc,
                           std::string type,
                           std::string extra);
  
  virtual void addAttValue(std::string name,
                             std::string value,
                             std::string showlabel);

  virtual void addAttValue(std::string name,
                             int         value,
                             std::string showlabel);

  virtual void addAttValue(std::string name,
                             float         value,
                             std::string showlabel);
  
  virtual void addPoint(float x, float y, float z);
  
  virtual void addInstance(std::string father, 
                             std::string type);

  virtual void addType(std::string father, 
                         std::string name, 
                         std::string desc, 
                         std::string info);

  void setTypeTree(HepRepTypeTree* tt){m_typeTree = tt; m_types.clear();};
  void setInstanceTree(HepRepInstanceTree* it){m_instanceTree = it; m_instances.clear();};

  HepRepType* getTypeBySequence(std::vector<int>& seq);
  HepRepInstance* getInstanceBySequence(std::vector<int>& seq);

  virtual void setSubinstancesNumber(std::string name, unsigned int n);

 private:
  typedef std::vector<int> sequenceCol;
  std::map<std::string,sequenceCol> m_types; 
  std::map<std::string,sequenceCol> m_instances; 

  bool isInType;

  sequenceCol m_actualType;
  //sequenceCol m_actualInstance;

  HepRepInstance* m_actualInstance;
  HepRepTypeTree* m_typeTree;
  HepRepInstanceTree* m_instanceTree;
};

#endif //CORBABUILDER_H

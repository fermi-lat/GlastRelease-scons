#ifndef HEPREP_IMPL_H
#define HEPREP_IMPL_H

#include <map>

#include "HepRepSvc/IRegistry.h"

#include "HepRep.hh"
#include "CorbaBuilder.h"

class HepRep_impl : virtual public POA_HepRep
{
 public:
  HepRep_impl();
  ~HepRep_impl();

  void addLayer (const char* layer);

  HepRepInstanceTree* getInstanceTreeTop(const char* instanceTreeName,
           const char* instanceTreeVersion);
  HepRepTypeTree* getTypeTree(const char* typeTreeName,
            const char* typeTreeVersion);
  HepRepInstanceTree* getInstances(const char* instanceTreeName,
           const char* instanceTreeVersion,
           const StringArray& typeNames);
  HepRepInstanceTree* getInstancesAfterAction(const char* instanceTreeName,
                const char* instanceTreeVersion,
                const StringArray& typeNames,
                const HepRepActionList& actions,
                CORBA::Boolean getPoints,
                CORBA::Boolean getDrawAtts,
                CORBA::Boolean getNonDrawAtts,
                const StringArray& invertAtt);
                
  StringArray* getLayerOrder();
  char* checkForException();

  void addInstanceByType(std::string instance, std::string type)
    {m_instanceTrees[instance] = type;};

  void clearInstances(){m_instanceTrees.clear();};

  void setRegistry(IRegistry* r);

private:
  std::map<std::string,HepRepTypeTree*> m_typeTrees;
  std::map<std::string,std::string> m_instanceTrees;

  IRegistry* m_registry;

  HepRepInstanceTree* _instanceTree;
  HepRepTypeTree* _typeTree;

  int _nTypes;
  const char** _typeNames;
  int _nLayers;
  const char** _layers;
  std::vector<std::string> m_layers;
  char* _lastException;

  CorbaBuilder* m_builder;
};
#endif

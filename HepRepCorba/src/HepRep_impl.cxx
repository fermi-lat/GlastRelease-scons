#include <iostream>

#include "HepRep_impl.hh"
#include "HepRepSvc/IFiller.h"


HepRep_impl::HepRep_impl ()
{
  m_builder = new CorbaBuilder();

  _nLayers = 0;
  addLayer("Geometry");
  addLayer("Event");
}


HepRep_impl::~HepRep_impl ()
{
  delete _typeTree;
  delete _instanceTree;

  delete _typeNames;
  delete _layers;
}

/*
void HepRep_impl::setHepRepSvc(IHepRepSvc* hrsv)
{
  m_hrsvc = hrsv;
  m_hrsvc->useBuilder(m_builder);

  typedef std::vector<IFiller*> fillerCol;
  std::vector<std::string>& types = m_hrsvc->getTypeTrees();  
  std::vector<std::string>::const_iterator i;  
  
  
  std::cout << "CORBA HepRep Object Ready " << std::endl;
  
  for(i=types.begin();i!=types.end();i++)
    {
      fillerCol::const_iterator j;
      HepRepTypeTree* typeTree = new HepRepTypeTree;
      typeTree->id.name = CORBA::string_dup((*i).c_str());
      typeTree->id.version = CORBA::string_dup("1.0");
      typeTree->types.length(0);
      
      m_builder->setTypeTree(typeTree);
      fillerCol temp = m_hrsvc->getFillersByType(*i);
      for(j=temp.begin(); j!=temp.end();j++)
        {
        (*j)->buildTypes();
        }
      
      m_typeTrees[(*i)]= typeTree;
    }
}
*/

HepRepInstanceTree* HepRep_impl::getInstanceTreeTop(const char* instanceTreeName,
						     const char* instanceTreeVersion)
{
  HepRepInstanceTree* instanceTreeTop = new HepRepInstanceTree;
  std::string type = m_registry->getTypeByInstance(instanceTreeName);  
  
  instanceTreeTop->id.name = CORBA::string_dup(instanceTreeName);
  instanceTreeTop->id.version = CORBA::string_dup(instanceTreeVersion);
  instanceTreeTop->typeTreeID.name = CORBA::string_dup(type.c_str());
  instanceTreeTop->typeTreeID.version = CORBA::string_dup("1.0");
  
  instanceTreeTop->instanceTreeIDs.length(0);
  instanceTreeTop->instances.length(0);


  for(unsigned int i=0;i < (m_registry->getDependencies(instanceTreeName)).size(); i++)
  {
    std::string temp = (m_registry->getDependencies(instanceTreeName))[i];
    int iType = instanceTreeTop->instanceTreeIDs.length()+1;
    instanceTreeTop->instanceTreeIDs.length(iType);
    instanceTreeTop->instanceTreeIDs[iType-1].name = CORBA::string_dup(temp.c_str());
    instanceTreeTop->instanceTreeIDs[iType-1].version = CORBA::string_dup("1.0");
  }
  
  return instanceTreeTop;
}

void HepRep_impl::setRegistry(IRegistry* r)
{
  m_registry = r;

  m_registry->useBuilder(m_builder);

  typedef std::vector<IFiller*> fillerCol;
  const std::vector<std::string>& types = m_registry->getTypeTrees();  
  std::vector<std::string>::const_iterator i;  
  
  
  std::cout << "CORBA HepRep Object Ready " << std::endl;
  
  for(i=types.begin();i!=types.end();i++)
    {
      fillerCol::const_iterator j;
      HepRepTypeTree* typeTree = new HepRepTypeTree;
      typeTree->id.name = CORBA::string_dup((*i).c_str());
      typeTree->id.version = CORBA::string_dup("1.0");
      typeTree->types.length(0);
           
      m_builder->setTypeTree(typeTree);
      fillerCol temp = m_registry->getFillersByType(*i);

      for(j=temp.begin(); j!=temp.end();j++)
        {
          (*j)->buildTypes();
        }
      
      m_typeTrees[(*i)]= typeTree;
    }
}

HepRepTypeTree* HepRep_impl::getTypeTree(const char* typeTreeName,
					 const char* typeTreeVersion)
{
  //  std::cout << "HepRep_impl::getTypeTree called for name " 
  //	    << typeTreeName << ", version " 
  //	    << typeTreeVersion << std::endl;

  // Trivial implementation since I am only supporting a single typeTree.

  HepRepTypeTree* typeTree = new HepRepTypeTree(*m_typeTrees[typeTreeName]);

  return typeTree;
}


HepRepInstanceTree* HepRep_impl::getInstances(const char* instanceTreeName,
					      const char* instanceTreeVersion,
					      const StringArray& typeNames)
{
  m_registry->useBuilder(m_builder);
  typedef std::vector<IFiller*> fillerCol;
  
  int nTypes = typeNames.length();

  HepRepInstanceTree* instanceTree = new HepRepInstanceTree;

  const fillerCol temp = 
    m_registry->getFillersByType(
				 m_registry->getTypeByInstance(instanceTreeName));  

  fillerCol::const_iterator j;  

  if (temp.size())
    {
      HepRepInstanceList tmp;
      tmp.length(0);
      instanceTree->id.name = CORBA::string_dup(instanceTreeName);
      instanceTree->id.version = CORBA::string_dup(instanceTreeVersion);
      instanceTree->typeTreeID.name = CORBA::string_dup((m_registry->getTypeByInstance(instanceTreeName)).c_str());
      instanceTree->typeTreeID.version = CORBA::string_dup("1.0");
      instanceTree->instances = tmp;
    }

  std::vector<std::string> names;
  for(unsigned int i=0;i<typeNames.length();i++)
  {
    names.push_back(CORBA::string_dup(typeNames[i]));
  }

  m_builder->setInstanceTree(instanceTree);  
  for(j=temp.begin(); j!=temp.end();j++)
    {
      (*j)->fillInstances(names);
    }

  return instanceTree;
}


HepRepInstanceTree* HepRep_impl::getInstancesAfterAction(const char* instanceTreeName,
                                                         const char* instanceTreeVersion,
                                                         const StringArray& typeNames,
                                                         const HepRepActionList& actions,
                                                         CORBA::Boolean getPoints,
                                                         CORBA::Boolean getDrawAtts,
                                                         CORBA::Boolean getNonDrawAtts,
                                                         const StringArray& invertAtt)
{
  return getInstances(instanceTreeName, instanceTreeVersion, typeNames);
}


StringArray* HepRep_impl::getLayerOrder()
{
  StringArray* layerOrder = new StringArray;
  layerOrder->length(m_layers.size());
    
  for (int iLayer=0; iLayer < m_layers.size(); iLayer++)
      (*layerOrder)[iLayer] = CORBA::string_dup(m_layers[iLayer].c_str());

  return layerOrder;
}


char* HepRep_impl::checkForException()
{
  return CORBA::string_dup(_lastException);
}


void HepRep_impl::addLayer(const char* layer)
{
// Called from HesHepRepReg to add layers in desired order.
//  _layers[_nLayers++] = layer;

  m_layers.push_back(layer);
}

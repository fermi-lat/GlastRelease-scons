#include <iostream>

#include "CorbaBuilder.h"


void CorbaBuilder::addAttDef(std::string name,
                             std::string desc,
                             std::string category,
                             std::string extra)
{
  HepRepType* type;

  type = getTypeBySequence(m_actualType);
  
  if (type)
    {
      int attNumber = type->attDefs.length() + 1;
      type->attDefs.length(attNumber);
      type->attDefs[attNumber-1].name     = CORBA::string_dup(name.c_str());
      type->attDefs[attNumber-1].desc     = CORBA::string_dup(desc.c_str());
      type->attDefs[attNumber-1].category = CORBA::string_dup(category.c_str());
      type->attDefs[attNumber-1].extra    = CORBA::string_dup(extra.c_str());
    }
}  

void CorbaBuilder::addAttValue(std::string name,
                               std::string value,
                               std::string showlabel)
{
  if (isInType)
    {
      HepRepType* type = getTypeBySequence(m_actualType);
      
      if(type)
        {
          int attNumber = type->attValues.length() + 1;
          type->attValues.length(attNumber);
          type->attValues[attNumber-1].name     = CORBA::string_dup(name.c_str());
          type->attValues[attNumber-1].value     <<= CORBA::string_dup(value.c_str());
          type->attValues[attNumber-1].showLabel     = 0;
        }
    }
  else
    {
      //HepRepInstance* instance = getInstanceBySequence(m_actualInstance);
      HepRepInstance* instance = m_actualInstance;      
      if(instance)
        {
          int attNumber2 = instance->attValues.length() + 1;
          instance->attValues.length(attNumber2);
          instance->attValues[attNumber2-1].name     = CORBA::string_dup(name.c_str());
          instance->attValues[attNumber2-1].value    <<= CORBA::string_dup(value.c_str());
          instance->attValues[attNumber2-1].showLabel     = 0;
        }
    } 
}

void CorbaBuilder::addAttValue(std::string name,
                               int value,
                               std::string showlabel)
{
  if (isInType)
    {
      HepRepType* type = getTypeBySequence(m_actualType);
      
      if(type)
        {
          int attNumber = type->attValues.length() + 1;
          type->attValues.length(attNumber);
          type->attValues[attNumber-1].name     = CORBA::string_dup(name.c_str());
          type->attValues[attNumber-1].value    <<= CORBA::ULong(value);
          type->attValues[attNumber-1].showLabel     = 0;
        }
    }
  else
    {
      //HepRepInstance* instance = getInstanceBySequence(m_actualInstance);
      HepRepInstance* instance = m_actualInstance;
      if(instance)
        {
          int attNumber2 = instance->attValues.length() + 1;
          instance->attValues.length(attNumber2);
          instance->attValues[attNumber2-1].name     = CORBA::string_dup(name.c_str());
          instance->attValues[attNumber2-1].value    <<= CORBA::ULong(value);
          instance->attValues[attNumber2-1].showLabel     = 0;
        }
    } 
}

void CorbaBuilder::addAttValue(std::string name,
                               float value,
                               std::string showlabel)
{
  if (isInType)
    {
      HepRepType* type = getTypeBySequence(m_actualType);
      
      if(type)
        {
          int attNumber = type->attValues.length() + 1;
          type->attValues.length(attNumber);
          type->attValues[attNumber-1].name     = CORBA::string_dup(name.c_str());
          type->attValues[attNumber-1].value      <<= CORBA::Double(value);
          type->attValues[attNumber-1].showLabel     = 0;
        }
    }
  else
    {
      //      HepRepInstance* instance = getInstanceBySequence(m_actualInstance);

      HepRepInstance* instance = m_actualInstance;
      
      if(instance)
        {
          int attNumber2 = instance->attValues.length() + 1;
          instance->attValues.length(attNumber2);
          instance->attValues[attNumber2-1].name     = CORBA::string_dup(name.c_str());
          instance->attValues[attNumber2-1].value    <<= CORBA::Double(value);
          instance->attValues[attNumber2-1].showLabel     = 0;
        }
    }   
}
  
void CorbaBuilder::addPoint(float x, float y, float z)
{
  //  HepRepInstance* instance;
   
  //  instance = getInstanceBySequence(m_actualInstance); 

  HepRepInstance* instance = m_actualInstance;

  int iType = instance->points.length()+1;
  instance->points.length(iType);
  instance->points[iType-1].x = x;
  instance->points[iType-1].y = y;
  instance->points[iType-1].z = z;
}
  
void CorbaBuilder::addInstance(std::string father, 
                               std::string type)
{
  /*
  unsigned int num;
  if (type == "ParticleCol")
    num = 3000;
  else
    num = 0;

  HepRepInstance::_tao_seq_HepRepInstance in(num);
  in.length(0);
  */

  if (father == "")
    {
      int iType = m_instanceTree->instances.length()+1;
      m_instanceTree->instances.length(iType);
      m_instanceTree->instances[iType-1].typeName = CORBA::string_dup(type.c_str());
      m_instanceTree->instances[iType-1].points.length(0);
      m_instanceTree->instances[iType-1].attValues.length(0);
      (m_instanceTree->instances[iType-1]).instances.length(0);
      
      std::vector<int> temp(1,iType-1);
      m_instances[type] = temp;
      //      m_actualInstance = temp;
      m_actualInstance = &m_instanceTree->instances[iType-1];
    }
  else
    {
      std::vector<int> temp2 = m_instances[father];

      HepRepInstance* instance;


      if(temp2.size())
        {
          instance = getInstanceBySequence(temp2);
        }

      if (instance)
        {
          int jType = instance->instances.length()+1;
          instance->instances.length(jType);


          instance->instances[jType-1].typeName = CORBA::string_dup(type.c_str());
          instance->instances[jType-1].points.length(0);
          instance->instances[jType-1].attValues.length(0);
          instance->instances[jType-1].instances.length(0);

          temp2.push_back(jType-1);
          m_instances[type] = temp2;
          //          m_actualInstance = temp2;

          m_actualInstance = &instance->instances[jType-1];
        }
    }  
  
  isInType = 0;
}

void CorbaBuilder::addType(std::string father, 
                           std::string name, 
                           std::string desc, 
                           std::string info)
{
    
  if (father == "")
    {
      int iType = m_typeTree->types.length()+1;
      m_typeTree->types.length(iType);
      m_typeTree->types[iType-1].name = CORBA::string_dup(name.c_str());
      m_typeTree->types[iType-1].desc = CORBA::string_dup(desc.c_str());
      m_typeTree->types[iType-1].infoURL = CORBA::string_dup(info.c_str());
      m_typeTree->types[iType-1].attDefs.length(0);
      m_typeTree->types[iType-1].attValues.length(0);
      m_typeTree->types[iType-1].types.length(0);

      
      std::vector<int> temp(1,iType-1);
      m_types[name] = temp;
      m_actualType = temp;
    }
  else
    {

      std::vector<int> temp2 = m_types[father];

      HepRepType* type;

      if(temp2.size())
        {
          type = getTypeBySequence(temp2);
        }

      
      if (type)
        {
          int jType = type->types.length()+1;
          type->types.length(jType);
          
          type->types[jType-1].name = CORBA::string_dup(name.c_str());
          type->types[jType-1].desc = CORBA::string_dup(desc.c_str());
          type->types[jType-1].infoURL = CORBA::string_dup(info.c_str());
          type->types[jType-1].attDefs.length(0);
          type->types[jType-1].attValues.length(0);
          type->types[jType-1].types.length(0);
                      
          temp2.push_back(jType-1);
          m_types[name] = temp2;
          m_actualType = temp2;
        }
    }

  isInType = 1;
}


HepRepType* CorbaBuilder::getTypeBySequence(std::vector<int>& seq)
{
  HepRepType* type = 0;
  
  if(seq.size())
    {
      type = &m_typeTree->types[seq[0]];          
      for(unsigned int i= 1; i<seq.size(); i++)
        {
          type = &(type->types[seq[i]]);              
        }
    }  

  return type;
}

HepRepInstance* CorbaBuilder::getInstanceBySequence(std::vector<int>& seq)
{
  HepRepInstance* instance = 0;
  
  if(seq.size())
    {
      instance = &m_instanceTree->instances[seq[0]];          
      for(unsigned int i= 1; i<seq.size(); i++)
        {
          instance = &(instance->instances[seq[i]]);              
        }
    }  

  return instance;
}

void CorbaBuilder::setSubinstancesNumber(std::string name, unsigned int i)
{
  std::vector<int> temp = m_instances[name];
  HepRepInstance* instance = getInstanceBySequence(temp);

  
  if (instance)
    {
//      HepRepInstanceList in(i);

//      instance->instances = in;
      
      instance->instances.replace(i,0,0,false);
    }
}


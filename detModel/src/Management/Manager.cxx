#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Management/Visitor.h"
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Management/ConstantsVisitor.h"
#include "detModel/Management/MaterialsVisitor.h"
#include "detModel/Management/Manager.h"
#include "detModel/Management/Builder.h"
#include "detModel/Gdd.h"
#include "detModel/Constants/Const.h"
#include "detModel/Constants/FloatConst.h"
#include "detModel/Constants/IntConst.h"
#include "detModel/Constants/DoubleConst.h"

namespace detModel{

  Manager* Manager::pointer = 0;

  Manager::~Manager(){
    delete manGdd;
    delete manBuilder;
  }

  Manager* Manager::getPointer()
  {
    if(pointer == 0)
      pointer = new Manager;
    return pointer;
  }

  void Manager::cleanGdd(){
    delete manGdd; 
    manGdd = new Gdd; 
    //  delete manBuilder;
  }

  void Manager::setBuilder(Builder* b)
  {
    manBuilder = b;
  
  }

  bool Manager::build(buildType b) 
  { 
    if (!manBuilder->parseFile(getNameFile())) return false;

    manGdd = manBuilder->getGdd();
  
    switch(b)
      {
      case all:
	manBuilder->buildSections(); 
	manBuilder->buildConstants(); 
	manBuilder->buildMaterials();
	manBuilder->buildIdDictionary(); 
	break;
      case constants:
	manBuilder->buildConstants(); 
	break;
      case sections:
	manBuilder->buildSections(); 
	break;
      case materials:
	manBuilder->buildMaterials();
	break;
      case idDictionary:
	manBuilder->buildIdDictionary();
	break;

      }
    return true;
  }

  void Manager::startVisitor(Visitor* v)
  {
    if (SectionsVisitor* sv = dynamic_cast<SectionsVisitor*>(v))
      {
      
	if (v->getRecursive())
	  manGdd->Accept(sv);
	else
	  manGdd->AcceptNotRec(sv);
      }
    else if(ConstantsVisitor* cv = dynamic_cast<ConstantsVisitor*>(v))
      {
	if (v->getRecursive())
	  manGdd->Accept(cv);
	else
	  manGdd->AcceptNotRec(cv);
      }
    else if(MaterialsVisitor* mv = dynamic_cast<MaterialsVisitor*>(v))
      {
	if (v->getRecursive())
	  manGdd->Accept(mv);
	else
	  manGdd->AcceptNotRec(mv);
      }
  }

  bool Manager::getNumericConstByName(std::string name, double* res)
  {
    Gdd* gdd = getGdd();
    Const* con = gdd->getConstByName(name);

    if (con) 
      {
	if (FloatConst* fl = dynamic_cast<FloatConst*>(con))
	  *res = fl->getValue();
	else if (DoubleConst* db = dynamic_cast<DoubleConst*>(con))
	  *res = db->getValue();
	else if (IntConst* in = dynamic_cast<IntConst*>(con))
	  *res = in->getValue();
	else return false;
	return true;
      }
    else return false;
  }

  bool Manager::getNumericConstByName(std::string name, int* res)
  {
    Gdd* gdd = getGdd();
    Const* con = gdd->getConstByName(name);

    if (con) 
      {
	if (IntConst* in = dynamic_cast<IntConst*>(con)) {
	  *res = in->getValue();
          return true;
        }
      }
    return false;
  }

}

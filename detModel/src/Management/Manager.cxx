#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Management/Visitor.h"
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Management/ConstantsVisitor.h"
#include "detModel/Management/Manager.h"
#include "detModel/Management/Builder.h"
#include "detModel/Gdd.h"

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

  void Manager::build(buildType b) 
  { 
    manBuilder->parseFile(getNameFile());
    manGdd = manBuilder->getGdd();
  
    switch(b)
      {
      case all:
	manBuilder->buildSections(); 
	manBuilder->buildConstants(); 
	break;
      case constants:
	manBuilder->buildConstants(); 
	break;
      case sections:
	manBuilder->buildSections(); 
	break;
      }
  }

  void Manager::startVisitor(Visitor* v)
  {
    /*
      switch(v->getType()){
      case SectionsVisitor:
    */
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
  }

}


#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Management/GDDvisitor.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Management/GDDconstantsVisitor.h"
#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDbuilder.h"

#include "detModel/GDD.h"


GDDmanager* GDDmanager::pointer = 0;

GDDmanager::~GDDmanager(){
  delete manGDD;
  delete manBuilder;
}

GDDmanager* GDDmanager::getPointer()
{
  if(pointer == 0)
    pointer = new GDDmanager;
  return pointer;
}

void GDDmanager::cleanGDD(){
  delete manGDD; 
  manGDD = new GDD; 
  //  delete manBuilder;
}

void GDDmanager::setBuilder(GDDbuilder* b)
{
  manBuilder = b;
  
}

void GDDmanager::build(buildType b) 
{ 
  manBuilder->parseFile(getNameFile());
  manGDD = manBuilder->getGDD();
  
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

void GDDmanager::startVisitor(GDDvisitor* v)
{
  /*
  switch(v->getType()){
  case GDDsectionsVisitor:
  */
  if (GDDsectionsVisitor* sv = dynamic_cast<GDDsectionsVisitor*>(v))
    {
      
      if (v->getRecursive())
	manGDD->Accept(sv);
      else
	manGDD->AcceptNotRec(sv);
    }
  else if(GDDconstantsVisitor* cv = dynamic_cast<GDDconstantsVisitor*>(v))
    {
      if (v->getRecursive())
	manGDD->Accept(cv);
      else
	manGDD->AcceptNotRec(cv);
    }
}




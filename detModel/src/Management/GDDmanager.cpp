#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Management/GDDvisitor.h"
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDbuilder.h"

#include "detModel/GDD.h"


GDDmanager* GDDmanager::pointer = 0;

GDDmanager::~GDDmanager(){
  delete manGDD;
}

GDDmanager* GDDmanager::getPointer()
{
  if(pointer == 0)
    pointer = new GDDmanager;
  return pointer;
}

void GDDmanager::setBuilder(GDDbuilder* b)
{
  manBuilder = b;
  manGDD = b->getGDD();
}

void GDDmanager::build(buildType b) 
{ 
  switch(b)
    {
    case all:
      manBuilder->buildConstants(); 
      manBuilder->buildSections(); 
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
  if (v->getType() == sectionsVisitor)
    {
      GDDsectionsVisitor* sv = static_cast<GDDsectionsVisitor*>(v);
      if (v->getRecursive())
	manGDD->Accept(sv);
      else
	manGDD->AcceptNotRec(sv);
    }
}



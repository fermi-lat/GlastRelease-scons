#include <string>
#include <iostream>
#include <fstream>
#include "vector.h"

#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDFakeBuilder.h"


#include "detModel/Sections/GDDsection.h"
#include "detModel/Sections/GDDbox.h"
#include "detModel/Sections/GDDcomposition.h"
#include "detModel/Sections/GDDposXYZ.h"
#include "detModel/Sections/GDDstack.h"
#include "detModel/Sections/GDDaxisPos.h"
#include "detModel/Sections/GDDaxisMPos.h"
#include "detModel/Sections/GDDchoice.h"
#include "detModel/Sections/GDDidField.h"
#include "detModel/Sections/GDDseg.h"

#include "detModel/Constants/GDDconstants.h"

#include "detModel/Constants/GDDintConst.h"
#include "detModel/Constants/GDDfloatConst.h"
#include "detModel/Constants/GDDdoubleConst.h"
#include "detModel/Constants/GDDstringConst.h"

GDDFakeBuilder::GDDFakeBuilder()
{

  GDDmanager* man = GDDmanager::getPointer();
  currentGDD = man->getGDD();

  currentGDD->setCVSid("fakeid");
  currentGDD->setDTDversion("fakeversion");

}



void GDDFakeBuilder::parseFile(char* nameFile){

  std::cout << "Fake parsing of the file" << std::endl; 

}  


///////////////////////////////////////////////////////////////////////////////////////////////
void GDDFakeBuilder::buildConstants(){

  std::cout << "Fake building of constants" << std::endl; 

  GDDconstants* ConstantsBranch = new GDDconstants();
  ConstantsBranch->setVersion("fakeMajor","fakeMinor");

  GDDconstCategory* cat= new GDDconstCategory;
  cat->setName("Fake Material List");
  cat->setOverview("I need all this stuff for testing purpouses ...");

  GDDstringConst* c1 = new GDDstringConst;
  c1->setName("BoxAMaterial");
  c1->setConstMeaning("mat");
  c1->setConstType(s);
  c1->setValue("Silicon");
  c1->setNote("I need this fake material");

  GDDstringConst* c2 = new GDDstringConst;
  c2->setName("BoxBMaterial");
  c2->setConstMeaning("mat");
  c2->setConstType(s);
  c2->setValue("W");
  c2->setNote("I need this fake material");

  cat->addConstant(c1);
  cat->addConstant(c2);

  ConstantsBranch->addConstantCategory(cat);

  currentGDD->setConstants(ConstantsBranch);
  currentGDD->buildConstantsMap();
}


void GDDFakeBuilder::buildSections()
{

  std::cout << "Fake building of constants" << std::endl; 

  // We build a new section
  GDDsection* s = new GDDsection("sec","1.0","17-07-01","RG","world");


  // We build some boxes to fill the section
  GDDbox* w = new GDDbox("world");

  w->setUnitLength("mm");
  w->setUnitAngle("deg");

  w->setX(500);
  w->setY(500);
  w->setZ(500);
  w->setMaterial("Vacuum");
  w->setSensitive(0);

  GDDbox* b1 = new GDDbox("boxA");

  b1->setUnitLength("mm");
  b1->setUnitAngle("deg");

  b1->setX(50);
  b1->setY(150);
  b1->setZ(100);
  b1->setMaterial("Silicon");
  b1->setSensitive(1);

  GDDbox* b2 = new GDDbox("boxB");

  b2->setUnitLength("mm");
  b2->setUnitAngle("deg");

  b2->setX(70);
  b2->setY(70);
  b2->setZ(70);
  b2->setMaterial("W");
  b2->setSensitive(0);


  // We build a composition to positionate this box
  GDDcomposition* c = new GDDcomposition();  

  c->setName("Comp");
  c->setEnvelope(w);

  GDDposXYZ* pos1 = new GDDposXYZ(posXYZ);

  pos1->setX(-50);
  pos1->setY(-50);
  pos1->setZ(0);

  pos1->setVolume(b1);

  GDDposXYZ* pos2 = new GDDposXYZ(posXYZ);

  pos2->setX(70);
  pos2->setY(40);
  pos2->setZ(20);

  pos2->setVolume(b2);

  c->addPosition(pos1);
  c->addPosition(pos2);

  s->setTopVolume(w);
  s->addVolume(w);
  s->addVolume(b1);
  s->addVolume(b2);
  s->addVolume(c);

  // We push the section in the currentGDD sections list
  currentGDD->getSections()->push_back(s);

  currentGDD->buildVolumeMap();

  currentGDD->buildBoundingBoxes();
}



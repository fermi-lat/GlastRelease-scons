#include <string>
#include <iostream>
#include <fstream>
#include <vector>

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

  GDDstringConst* c3 = new GDDstringConst;
  c3->setName("BoxCMaterial");
  c3->setConstMeaning("mat");
  c3->setConstType(s);
  c3->setValue("CsI");
  c3->setNote("I need this fake material");

  GDDstringConst* c4 = new GDDstringConst;
  c4->setName("BoxDMaterial");
  c4->setConstMeaning("mat");
  c4->setConstType(s);
  c4->setValue("C");
  c4->setNote("I need this fake material");

  cat->addConstant(c1);
  cat->addConstant(c2);
  cat->addConstant(c3);
  cat->addConstant(c4);

  ConstantsBranch->addConstantCategory(cat);

  currentGDD->setConstants(ConstantsBranch);
  currentGDD->buildConstantsMap();
}


void GDDFakeBuilder::buildSections()
{

  std::cout << "Fake building of sections" << std::endl; 

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

  b1->setX(40);
  b1->setY(40);
  b1->setZ(40);
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

  GDDbox* b3 = new GDDbox("boxC");

  b3->setUnitLength("mm");
  b3->setUnitAngle("deg");

  b3->setX(70);
  b3->setY(70);
  b3->setZ(70);
  b3->setMaterial("CsI");
  b3->setSensitive(0);

  GDDbox* b4 = new GDDbox("boxD");

  b4->setUnitLength("mm");
  b4->setUnitAngle("deg");

  b4->setX(70);
  b4->setY(70);
  b4->setZ(70);
  b4->setMaterial("C");
  b4->setSensitive(0);


  // We build a composition to positionate this box
  GDDcomposition* c = new GDDcomposition();  

  c->setName("Comp");
  c->setEnvelope(w);

  GDDposXYZ* pos1 = new GDDposXYZ(posXYZ);

  pos1->setX(0);
  pos1->setY(0);
  pos1->setZ(0);

  pos1->setVolume(b1);

  GDDposXYZ* pos2 = new GDDposXYZ(posXYZ);

  pos2->setX(100);
  pos2->setY(0);
  pos2->setZ(0);

  pos2->setVolume(b2);

  GDDposXYZ* pos3 = new GDDposXYZ(posXYZ);

  pos3->setX(0);
  pos3->setY(100);
  pos3->setZ(0);

  pos3->setVolume(b3);

  GDDposXYZ* pos4 = new GDDposXYZ(posXYZ);

  pos4->setX(0);
  pos4->setY(0);
  pos4->setZ(100);

  pos4->setVolume(b4);

  c->addPosition(pos1);
  c->addPosition(pos2);
  c->addPosition(pos3);
  c->addPosition(pos4);

  s->setTopVolume(w);
  s->addVolume(w);
  s->addVolume(b1);
  s->addVolume(b2);
  s->addVolume(b3);
  s->addVolume(b4);
  s->addVolume(c);

  // We push the section in the currentGDD sections list
  currentGDD->getSections()->push_back(s);

  currentGDD->buildVolumeMap();

  currentGDD->buildBoundingBoxes();
}



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
  c1->setValue("Silicon");
  c1->setNote("I need this fake material");

  GDDstringConst* c2 = new GDDstringConst;
  c2->setName("BoxBMaterial");
  c2->setConstMeaning("mat");
  c2->setValue("W");
  c2->setNote("I need this fake material");

  GDDstringConst* c3 = new GDDstringConst;
  c3->setName("BoxCMaterial");
  c3->setConstMeaning("mat");
  c3->setValue("CsI");
  c3->setNote("I need this fake material");

  GDDstringConst* c4 = new GDDstringConst;
  c4->setName("BoxDMaterial");
  c4->setConstMeaning("mat");
  c4->setValue("C");
  c4->setNote("I need this fake material");

  GDDstringConst* c5 = new GDDstringConst;
  c5->setName("WorldMaterial");
  c5->setConstMeaning("mat");
  c5->setValue("Vacuum");
  c5->setNote("I need this fake material");

  cat->addConstant(c1);
  cat->addConstant(c2);
  cat->addConstant(c3);
  cat->addConstant(c4);
  cat->addConstant(c5);

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

  w->getUnits()->setUnitLength("mm");
  w->getUnits()->setUnitAngle("deg");

  w->setX(500);
  w->setY(500);
  w->setZ(500);
  w->setMaterial("Vacuum");
  w->setSensitive(0);

  GDDbox* b1 = new GDDbox("boxA");

  b1->getUnits()->setUnitLength("mm");
  b1->getUnits()->setUnitAngle("deg");

  b1->setX(40);
  b1->setY(40);
  b1->setZ(40);
  b1->setMaterial("Silicon");
  b1->setSensitive(1);

  GDDbox* b2 = new GDDbox("boxB");

  b2->getUnits()->setUnitLength("mm");
  b2->getUnits()->setUnitAngle("deg");

  b2->setX(70);
  b2->setY(70);
  b2->setZ(70);
  b2->setMaterial("W");
  b2->setSensitive(0);

  GDDbox* b3 = new GDDbox("boxC");

  b3->getUnits()->setUnitLength("mm");
  b3->getUnits()->setUnitAngle("deg");

  b3->setX(70);
  b3->setY(70);
  b3->setZ(70);
  b3->setMaterial("CsI");
  b3->setSensitive(0);

  GDDbox* b4 = new GDDbox("boxD");

  b4->getUnits()->setUnitLength("mm");
  b4->getUnits()->setUnitAngle("deg");

  b4->setX(70);
  b4->setY(70);
  b4->setZ(70);
  b4->setMaterial("C");
  b4->setSensitive(0);


  // We build a composition to positionate some boxes
  GDDcomposition* c = new GDDcomposition();    
  c->setName("Comp");
  c->setEnvelope(w);

  // We build a stack to positionate some boxes
  GDDstack* st = new GDDstack(GDDstack::zDir);  
  st->setName("Stack");

  GDDstack* st2 = new GDDstack(GDDstack::xDir);  
  st2->setName("Stack2");

  GDDaxisMPos* ampos1 = new GDDaxisMPos(GDDaxisPos::zDir);
  ampos1->setVolume(b1);
  ampos1->setNcopy(3);
  ampos1->setGap(20);

  st->addPosition(ampos1);

  GDDaxisMPos* ampos2 = new GDDaxisMPos(GDDaxisPos::xDir);
  ampos2->setVolume(st);
  ampos2->setNcopy(3);
  ampos2->setGap(40);
  st2->addPosition(ampos2);  


  GDDposXYZ* pos1 = new GDDposXYZ(0,0,0);
  pos1->setVolume(st2);
  
  c->addPosition(pos1);

  s->setTopVolume(c);
  s->setTopVolumeRef(w->getName());
  s->addVolume(w);
  s->addVolume(b1);
  s->addVolume(b2);
  s->addVolume(b3);
  s->addVolume(b4);
  s->addVolume(c);
  s->addVolume(st);
  s->addVolume(st2);

  // We push the section in the currentGDD sections list
  currentGDD->addSection(s);

  currentGDD->buildVolumeMap();

  currentGDD->buildBoundingBoxes();

}



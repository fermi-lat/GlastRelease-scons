#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "detModel/Management/Manager.h"
#include "detModel/Management/FakeBuilder.h"
#include "detModel/Sections/Section.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Composition.h"
#include "detModel/Sections/PosXYZ.h"
#include "detModel/Sections/Stack.h"
#include "detModel/Sections/AxisMPos.h"
#include "detModel/Constants/Constants.h"
#include "detModel/Constants/StringConst.h"

namespace detModel{

  FakeBuilder::FakeBuilder()
  {

    Manager* man = Manager::getPointer();
    currentGdd = man->getGdd();

    currentGdd->setCVSid("fakeid");
    currentGdd->setDTDversion("fakeversion");

  }



  void FakeBuilder::parseFile(char* nameFile){
  
    std::cout << "Fake parsing of the file " << nameFile << std::endl; 

  }  


  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FakeBuilder::buildConstants(){

    std::cout << "Fake building of constants" << std::endl; 
  
    Constants* ConstantsBranch = new Constants();
    ConstantsBranch->setVersion("fakeMajor","fakeMinor");

    ConstCategory* cat= new ConstCategory;
    cat->setName("Fake Material List");
    cat->setOverview("I need all this stuff for testing purpouses ...");

    StringConst* c1 = new StringConst;
    c1->setName("BoxAMaterial");
    c1->setConstMeaning("mat");
    c1->setValue("Silicon");
    c1->setNote("I need this fake material");

    StringConst* c2 = new StringConst;
    c2->setName("BoxBMaterial");
    c2->setConstMeaning("mat");
    c2->setValue("W");
    c2->setNote("I need this fake material");

    StringConst* c3 = new StringConst;
    c3->setName("BoxCMaterial");
    c3->setConstMeaning("mat");
    c3->setValue("CsI");
    c3->setNote("I need this fake material");

    StringConst* c4 = new StringConst;
    c4->setName("BoxDMaterial");
    c4->setConstMeaning("mat");
    c4->setValue("C");
    c4->setNote("I need this fake material");

    StringConst* c5 = new StringConst;
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

    currentGdd->setConstants(ConstantsBranch);
    currentGdd->buildConstantsMap();
  }


  void FakeBuilder::buildSections()
  {

    std::cout << "Fake building of sections" << std::endl; 

    // We build a new section
    Section* s = new Section("sec","1.0","17-07-01","RG","","world");
  
    // We build some boxes to fill the section
    Box* w = new Box("world");

    w->getUnits()->setUnitLength("mm");
    w->getUnits()->setUnitAngle("deg");

    w->setX(500);
    w->setY(500);
    w->setZ(500);
    w->setMaterial("Vacuum");
    w->setSensitive(0);

    Box* b1 = new Box("boxA");

    b1->getUnits()->setUnitLength("mm");
    b1->getUnits()->setUnitAngle("deg");

    b1->setX(40);
    b1->setY(40);
    b1->setZ(40);
    b1->setMaterial("Silicon");
    b1->setSensitive(1);

    Box* b2 = new Box("boxB");

    b2->getUnits()->setUnitLength("mm");
    b2->getUnits()->setUnitAngle("deg");

    b2->setX(70);
    b2->setY(70);
    b2->setZ(70);
    b2->setMaterial("W");
    b2->setSensitive(0);

    Box* b3 = new Box("boxC");

    b3->getUnits()->setUnitLength("mm");
    b3->getUnits()->setUnitAngle("deg");

    b3->setX(70);
    b3->setY(70);
    b3->setZ(70);
    b3->setMaterial("CsI");
    b3->setSensitive(0);

    Box* b4 = new Box("boxD");

    b4->getUnits()->setUnitLength("mm");
    b4->getUnits()->setUnitAngle("deg");

    b4->setX(70);
    b4->setY(70);
    b4->setZ(70);
    b4->setMaterial("C");
    b4->setSensitive(0);


    // We build a composition to positionate some boxes
    Composition* c = new Composition();    
    c->setName("Comp");
    c->setEnvelope(w);

    // We build a stack to positionate some boxes
    Stack* st = new Stack(Stack::zDir);  
    st->setName("Stack");

    Stack* st2 = new Stack(Stack::xDir);  
    st2->setName("Stack2");

    AxisMPos* ampos1 = new AxisMPos(Stack::zDir);
    ampos1->setVolume(b1);
    ampos1->setNcopy(3);
    ampos1->setGap(20);

    st->addPosition(ampos1);

    AxisMPos* ampos2 = new AxisMPos(Stack::xDir);
    ampos2->setVolume(st);
    ampos2->setNcopy(3);
    ampos2->setGap(40);
    st2->addPosition(ampos2);  


    PosXYZ* pos1 = new PosXYZ(0,0,0);
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

    // We push the section in the currentGdd sections list
    currentGdd->addSection(s);

    currentGdd->buildVolumeMap();

    currentGdd->buildBoundingBoxes();

  }

}

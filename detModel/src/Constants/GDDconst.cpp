#include "detModel/Constants/GDDconst.h"
#include "detModel/Management/GDDconstantsVisitor.h"


void GDDconst::Accept(GDDconstantsVisitor* v){
  v->visitConst(this);
}

void GDDconst::AcceptNotRec(GDDconstantsVisitor* v){
  v->visitConst(this);
}


void GDDconst::setConstMeaning(std::string puType){
  if (puType=="length")uType=length;
  else if (puType=="angle")uType=angle;
  else if (puType=="count")uType=coun;
  else if (puType=="mat")uType=mat;
  else if (puType=="none")uType=none;
}

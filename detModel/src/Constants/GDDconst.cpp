#include "detModel/Constants/GDDconst.h"


void GDDconst::setConstMeaning(string puType){
  if (puType=="length")uType=length;
  else if (puType=="angle")uType=angle;
  else if (puType=="count")uType=coun;
  else if (puType=="mat")uType=mat;
  else if (puType=="none")uType=none;
}

#include "detModel/Constants/Const.h"
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{

  void Const::setConstMeaning(std::string puType){
    if (puType=="length")uType=length;
    else if (puType=="angle")uType=angle;
    else if (puType=="count")uType=coun;
    else if (puType=="mat")uType=mat;
    else if (puType=="energy")uType=energy;
    else if (puType=="none")uType=none;
  }

}

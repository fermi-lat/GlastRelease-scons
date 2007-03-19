//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Trigger window parameters in GEM configuration
//
// Environment:
//      Software developed for GLAST.
//
// Author List:
//      Martin Kocian
//
// Copyright Information:
//      Copyright (C) 2005      Stanford Linear Accelerator Center
//
//---------------------------------------------------------------------------

#include "configData/gem/TrgWindowParams.h"

TrgWindowParams::TrgWindowParams(){
  _winwidth=0;
  _windowmask=0;
}

int TrgWindowParams::windowMask()const{
  return (int)_windowmask;
}

int TrgWindowParams::windowWidth()const{
  return _winwidth;
}

bool TrgWindowParams::maskIncludes(TrgConditions::conditions co){
  return (bool)((1<<co)&_windowmask);
}

void TrgWindowParams::setWindowMask(int m){
  _windowmask=(char)m;
}

void TrgWindowParams::setWindowWidth(int w){
  _winwidth=w&0x1f;
}

void TrgWindowParams::addToOpenMask(TrgConditions::conditions cond){
  _windowmask|=(1<<cond);
}

std::ostream& operator <<(std::ostream &os, const TrgWindowParams& tl){
  os<<"Window open mask: ";
  for (int i=0;i<NCOND;i++){
    if ((1<<i)&tl.windowMask())os<<TrgConditions::condName((TrgConditions::conditions)i)<<" ";
  }
  os<<std::endl;
  os<<"Window width: "<<tl.windowWidth()<<std::endl;
}

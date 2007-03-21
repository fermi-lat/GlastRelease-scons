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
#include "TTree.h"

TrgWindowParams::TrgWindowParams()
  :ConfigBranch("window",'i',ChannelKey(1)){
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

// Attach this value to a TTree
void TrgWindowParams::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();

  std::string branchNameWidth = branchName; branchNameWidth += "_width";
  std::string leafNameWidth = branchNameWidth;
  setLeafSuffix(leafNameWidth);
  tree.Branch(branchNameWidth.c_str(),(void*)(&_winwidth),leafNameWidth.c_str());
  
  std::string branchNameMask = branchName; branchNameMask += "_open_mask";
  std::string leafNameMask = branchNameMask;
  leafNameMask += "/B";
  tree.Branch(branchNameMask.c_str(),(void*)(&_windowmask),leafNameMask.c_str());
}

void TrgWindowParams::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();

  std::string branchNameWidth = branchName; branchNameWidth += "_width";
  tree.SetBranchAddress(branchNameWidth.c_str(),(void*)(&_winwidth));
  
  std::string branchNameMask = branchName; branchNameMask += "_open_mask";
  tree.SetBranchAddress(branchNameMask.c_str(),(void*)(&_windowmask));
}


std::ostream& operator <<(std::ostream &os, const TrgWindowParams& tl){
  os<<"Window open mask: ";
  for (int i=0;i<NCOND;i++){
    if ((1<<i)&tl.windowMask())os<<TrgConditions::condName((TrgConditions::conditions)i)<<" ";
  }
  os<<std::endl;
  os<<"Window width: "<<tl.windowWidth()<<std::endl;
}

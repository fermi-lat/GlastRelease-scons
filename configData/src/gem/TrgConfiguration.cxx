//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      A GEM configuration register class 
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

#include "configData/gem/TrgConfiguration.h"
#include <iomanip>
#include "TTree.h"

TrgConfiguration::TrgConfiguration():
  ConfigBranch("configuration",'i',ChannelKey(1)),
  _configuration(0){}

TrgConfiguration::TrgConfiguration(unsigned long c):
  ConfigBranch("configuration",'i',ChannelKey(1)),
  _configuration(c){}

void TrgConfiguration::setConfiguration(unsigned long c){
  _configuration=c&0x7f;
}

void TrgConfiguration::useAcdAsTrigger(bool on){
  if (on)_configuration|=1;
  else _configuration&=0x7e;
}
bool TrgConfiguration::acdUsedAsTrigger() const{
  return _configuration&0x1;
}

// Attach this value to a TTree
void TrgConfiguration::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  std::string leafName = branchName;
  setLeafSuffix(leafName);
  tree.Branch(branchName.c_str(),(void*)(&_configuration),leafName.c_str());
}

void TrgConfiguration::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  tree.SetBranchAddress(branchName.c_str(),(void*)(&_configuration));
}

std::ostream& operator <<(std::ostream& os,  const TrgConfiguration& tc){
  if(tc.acdUsedAsTrigger())os<<"ACD used as trigger"<<std::endl;
  else os<<"ACD used as veto"<<std::endl;
  os<<"Configuration register: "<<std::hex<<tc.configuration()<<std::dec<<std::endl;
  return os;
}

//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Periodic Trigger Parameters in GEM configuration
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

#include "configData/gem/TrgPeriodicTrigger.h"

#include "TTree.h"

TrgPeriodicTrigger::TrgPeriodicTrigger()
  :ConfigBranch("periodic",'i',ChannelKey(1)){
  _prescale=0;
  _limit=0;
  _free=false;
  _onepps=false;
}

int TrgPeriodicTrigger::prescale() const{
  return _prescale;
}
int TrgPeriodicTrigger::limit() const{
  return _limit;
}
bool TrgPeriodicTrigger::onePPS() const{
  return _onepps;
}
bool TrgPeriodicTrigger::freeRunning() const{
  return _free;
}
void TrgPeriodicTrigger::setPrescale(int p){
  _prescale=p&0xffffff;
}
void TrgPeriodicTrigger::setLimit(int l){
  _limit=l&0xffff;
}
void TrgPeriodicTrigger::setOnePPS(bool o){
  _onepps=o;
}
void TrgPeriodicTrigger::setFreeRunning(bool o){
  _free=o;
}

// Reset the cached and output values
void TrgPeriodicTrigger::reset() {
  _prescale = 0;
  _limit = 0;
  _free = false;
  _onepps = false;
}

// Attach this value to a TTree
void TrgPeriodicTrigger::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();

  std::string branchNamePres = branchName; branchNamePres += "_prescale";
  std::string leafNamePres = branchNamePres;
  setLeafSuffix(leafNamePres);
  tree.Branch(branchNamePres.c_str(),(void*)(&_prescale),leafNamePres.c_str());
  
  std::string branchNameLimit = branchName; branchNameLimit += "_limit";
  std::string leafNameLimit = branchNameLimit;
  setLeafSuffix(leafNameLimit);
  tree.Branch(branchNameLimit.c_str(),(void*)(&_limit),leafNameLimit.c_str());

  std::string branchNameFree = branchName; branchNameFree += "_free";
  std::string leafNameFree = branchNameFree;
  leafNameFree += "/B";
  tree.Branch(branchNameFree.c_str(),(void*)(&_free),leafNameFree.c_str());

  std::string branchNameOnepps = branchName; branchNameOnepps += "_onepps";
  std::string leafNameOnepps = branchNameOnepps;
  leafNameOnepps += "/B";
  tree.Branch(branchNameOnepps.c_str(),(void*)(&_onepps),leafNameOnepps.c_str());
}

void TrgPeriodicTrigger::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();

  std::string branchNamePres = branchName; branchNamePres += "_prescale";
  tree.SetBranchAddress(branchNamePres.c_str(),(void*)(&_prescale));
  
  std::string branchNameLimit = branchName; branchNameLimit += "_limit";
  tree.SetBranchAddress(branchNameLimit.c_str(),(void*)(&_limit));

  std::string branchNameFree = branchName; branchNameFree += "_free";
  tree.SetBranchAddress(branchNameFree.c_str(),(void*)(&_free));

  std::string branchNameOnepps = branchName; branchNameOnepps += "_onepps";
  tree.SetBranchAddress(branchNameOnepps.c_str(),(void*)(&_onepps));
}

std::ostream& operator <<(std::ostream &os, const TrgPeriodicTrigger& tl){

  os<<"Periodic prescale: "<<tl.prescale()<<std::endl;
  os<<"1PPS: ";
  if (tl.onePPS())os<<"True"<<std::endl;
  else os<<"False"<<std::endl;
  os<<"Free running: ";
  if (tl.freeRunning())os<<"True"<<std::endl;
  else os<<"False"<<std::endl;
  os<<"Periodic limit: "<<tl.limit()<<std::endl;
  return os;
}

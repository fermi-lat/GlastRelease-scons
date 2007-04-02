//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      A trigger engine configuration class
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


#include "configData/gem/TrgEngine.h"
#include "TTree.h"

TrgEngine::TrgEngine()
  :ConfigBranch("engine",'i',ChannelKey(16)){
  clear();
}

void TrgEngine::setEngine(unsigned idx, unsigned long en){
  _engine[idx] =en;
}
unsigned long TrgEngine::engine(unsigned idx) const {
  return _engine[idx];
}
void TrgEngine::setCalstrobe(unsigned idx,bool on){
  if(on)_engine[idx]|=1<<16;
  else _engine[idx]&=~(1<<16);
}
void TrgEngine::setTack(unsigned idx,bool on){
  if(on)_engine[idx]|=1<<19;
  else _engine[idx]&=~(1<<19);
}
void TrgEngine::setZeroSuppression(unsigned idx,bool on){
  if(on)_engine[idx]|=1<<21;
  else _engine[idx]&=~(1<<21);
}
void TrgEngine::setFourRangeReadout(unsigned idx,bool on){
  if(on)_engine[idx]|=1<<20;
  else _engine[idx]&=~(1<<20);

}
void TrgEngine::setInhibit(unsigned idx,bool on){
  if(on)_engine[idx]|=1<<15;
  else _engine[idx]&=~(1<<15);
}
void TrgEngine::setDestination(unsigned idx,int dest){
  _engine[idx]&=~(0x3f<<25);
  _engine[idx]|=((dest&0x3f)<<25);
}
void TrgEngine::setMarker(unsigned idx,int marker){
  _engine[idx]&=~(0x7<<22);
  _engine[idx]|=((marker&0x7)<<22);
}
void TrgEngine::setPrescale(unsigned idx,int num){
  _engine[idx]&=~0xff;
  _engine[idx]|=num&0xff;
}
bool TrgEngine::calstrobe(unsigned idx) const{
  return _engine[idx]&(1<<16)? 1: 0;
}
bool TrgEngine::tack(unsigned idx) const{
  return _engine[idx]&(1<<19)? 1: 0;
}
bool TrgEngine::inhibited(unsigned idx) const{
  return _engine[idx]&(1<<15)? 1: 0;
}
bool TrgEngine::zeroSuppression(unsigned idx) const{
  return _engine[idx]&(1<<21)? 1: 0;
}
bool TrgEngine::fourRangeReadout(unsigned idx) const{
  return _engine[idx]&(1<<20)? 1: 0;
}
int TrgEngine::destination(unsigned idx) const{
  return (_engine[idx]>>25)&0x3f;
}
int TrgEngine::marker(unsigned idx) const{
  return (_engine[idx]>>22)&0x7;
}
int TrgEngine::prescale(unsigned idx) const{
  return _engine[idx]&0xff;
}
void TrgEngine::clear(){
  for (int i(0); i < 16; i++ ) {
    _engine[i]=0;
  }
} 

// Attach this value to a TTree
void TrgEngine::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  std::string leafName = branchName;
  setLeafSuffix(leafName);
  tree.Branch(branchName.c_str(),(void*)(_engine),leafName.c_str());
}

void TrgEngine::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  tree.SetBranchAddress(branchName.c_str(),(void*)(_engine));
}

std::ostream& operator <<(std::ostream& os, const TrgEngine& tc){
  for (int i(0); i < 16; i++ ) {
    os<<"Engine "<<i<<": ";
    if (tc.calstrobe(i))os<<"Calstrobe ";
    if (tc.tack(i))os<<"Tack ";
    if (tc.inhibited(i))os<<"Inhibit ";
    os<<"; ";
    if (tc.prescale(i)!=0)os<<"Prescale="<<tc.prescale(i)<<" ";
    if (tc.marker(i)!=0)os<<"Marker="<<tc.marker(i)<<" ";
    if (tc.destination(i)!=0)os<<"Destination="<<tc.destination(i)<<" ";
    if (tc.zeroSuppression(i))os<<"0-suppression ";
    if (tc.fourRangeReadout(i))os<<"4-range-readout ";
    os<<std::endl;
  }
  return os;
}
    


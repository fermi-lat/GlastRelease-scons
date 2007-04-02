//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      GEM LUT class
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

#include "configData/gem/TrgConfLUT.h"
#include <iomanip>
#include "TTree.h"

TrgConfLUT::TrgConfLUT()
  :ConfigBranch("conditions",'i',ChannelKey(NREG)){
    for (int i=0;i<NREG;i++) _lut[i]=0;
  }  

int TrgConfLUT::engineNumber(int i) const{
  if (i<0 || i>NENTRIES)return -1;
  unsigned long ro;
  ro=_lut[i/CONDSREG];
  ro=ro>>((i%CONDSREG)*NBITS);
  ro&=CONDMASK;  
  return ro;
}

void TrgConfLUT::setEngineNumber(int i, int ent){
  unsigned long mask;
  mask=CONDMASK<<((i%CONDSREG)*NBITS);
  mask=~mask;
  _lut[i/CONDSREG]&=mask;
  _lut[i/CONDSREG]|=(ent&CONDMASK)<<((i%CONDSREG)*NBITS);
}

std::vector<unsigned long> TrgConfLUT::conditionsInEngine(int num)const{
  std::vector<unsigned long> retvector;
  if (num<0 || num>15)return retvector;
  for (unsigned i=0;i<NENTRIES;i++){
    if(engineNumber(i)==num)retvector.push_back(i);
  }
  return retvector;

}

void TrgConfLUT::reset() {
  for ( int i(0); i < NREG; i++ ) {
    _lut[i] = 0;
  }
}

// Attach this value to a TTree
void TrgConfLUT::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  std::string leafName = branchName;
  setLeafSuffix(leafName);
  tree.Branch(branchName.c_str(),(void*)(_lut),leafName.c_str());
}

void TrgConfLUT::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  tree.SetBranchAddress(branchName.c_str(),(void*)(_lut));
}


std::ostream& operator <<(std::ostream& os, const TrgConfLUT& tc){
  os<<"Lookup-table registers:"<<std::endl;
  os<<std::hex<<std::setfill('0');
  for (int i=0;i<NREG;i++){
    os<<std::setw(CONDSREG)<<tc.reg(i)<<" ";
    if (i%CONDSREG==CONDSREG - 1)os<<std::endl;
  }
  os<<std::dec;
  return os;
}

//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      A GEM ROI register class 
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

#include "configData/gem/TrgRoi.h"
#include <iomanip>
#include "commonRootData/idents/AcdId.h"
#include "TTree.h"
#include <algorithm> 

TrgRoi::TrgRoi()
  :ConfigBranch("r",'i',ChannelKey(54)){  
  clear();
}

void TrgRoi::clear(){
  for (int i=0;i<54;i++){
    _roi[i]=0;
  }
}

void TrgRoi::setRoiRegister(int i, unsigned long val){
  if (i<0 || i>53)return;
  _roi[i]=val;
}

unsigned long TrgRoi::roiRegister(int i) const{
  if (i<0 || i>53)return 0;
  return _roi[i];
}

std::vector<unsigned long> TrgRoi::roiFromIndex(int index) const{
  std::vector<unsigned long> retvec;
  if(index>107 || index<0)return retvec;
  int reg=index/2;
  unsigned bitshift=(index%2)*16;
  for (unsigned i=0;i<16;i++){
    if (_roi[reg]&(1<<(bitshift+i)))retvec.push_back(i);
  }
  return retvec;
}

std::vector<unsigned long> TrgRoi::roiFromName(int name) const{
  int index=AcdId::indexFromTile(name);
  return roiFromIndex(index);
}

std::vector<unsigned long> TrgRoi::tilesByIndexInRoi(int roi) const{
  std::vector<unsigned long> retvec; 
  if (roi<0 || roi> 15)return retvec;
  for (int i=0;i<54;i++){
    if(_roi[i]&(1<<roi))retvec.push_back(i*2);
    if(_roi[i]&(1<<(roi+16)))retvec.push_back(i*2+1);
  } 
  return retvec;
}

std::vector<unsigned long> TrgRoi::tilesByNameInRoi(int roi) const{
  std::vector<unsigned long>retvec=tilesByIndexInRoi(roi);
  std::transform(retvec.begin(),retvec.end(),retvec.begin(),&AcdId::tileFromIndex);
  return retvec;
}

bool TrgRoi::isTileByIndexInRoi(unsigned tile, int roi) const{
  if (roi>15 || roi<0) return false;
  std::vector<unsigned long>roivec=tilesByIndexInRoi(roi);
  if (std::find(roivec.begin(),roivec.end(),tile)!=roivec.end())return true;
  else return false;
}

bool TrgRoi::isTileByNameInRoi(unsigned tile, int roi) const{
  if (roi>15 || roi<0) return false;
  std::vector<unsigned long>roivec=tilesByNameInRoi(roi);
  if (std::find(roivec.begin(),roivec.end(),tile)!=roivec.end())return true;
  else return false;
}
  

// Attach this value to a TTree
void TrgRoi::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  std::string leafName = branchName;
  setLeafSuffix(leafName);
  tree.Branch(branchName.c_str(),(void*)(_roi),leafName.c_str());
}

void TrgRoi::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  tree.SetBranchAddress(branchName.c_str(),(void*)(_roi));
}

std::ostream& operator <<(std::ostream& os,  const TrgRoi& tc){
  os<<"ROI: Tiles"<<std::endl;
  for (int i=0;i<16;i++){
    os<<std::setw(2)<<i<<":  ";
    std::vector<unsigned long> roivec=tc.tilesByNameInRoi(i);
    for (unsigned j=0;j<roivec.size();j++){
      os<<std::setfill('0')<<std::setw(3)<<roivec[j]<<" ";
    }
    os<<std::endl;
  }
  return os;
}


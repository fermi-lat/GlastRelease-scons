//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Disabled channels in GEM configuration
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

#include "configData/gem/TrgDisabledChannels.h"
#include "commonRootData/idents/AcdId.h"
#include <iomanip>
#include "TTree.h"

TrgDisabledChannels::TrgDisabledChannels()
  :ConfigBranch("",'i',ChannelKey(1)){
  clear();
}
void TrgDisabledChannels::clear(){
  for (int i=0;i<12;i++){
    _tiles[i]=0;
  }
  _towers[0]=_towers[1]=_towers[2]=_towers[3]=0;
  _cno=0;
  _busy=0;
  _external=0;


}

void TrgDisabledChannels::setTowerRegister(unsigned i, unsigned val){
  if (i<4)_towers[i]=val&0xfff;
}
unsigned TrgDisabledChannels::towerRegister(unsigned i){
  if (i<4)return _towers[i];
  else return 0;
}
void TrgDisabledChannels::setTileRegister(unsigned i, unsigned val){
  if (i<12)_tiles[i]=val&0x3ffff;
}
unsigned TrgDisabledChannels::tileRegister(unsigned i){
  if (i<12)return _tiles[i];
  else return 0;
}
void TrgDisabledChannels::setCnoRegister( unsigned val){
  _cno=val&0xfff;
}
void TrgDisabledChannels::setBusyRegister( unsigned val){
  _busy=val&0xffff;
}
void TrgDisabledChannels::setExternal( unsigned val){
  _external=(bool)val;
}

void TrgDisabledChannels::setTowerMask(unsigned tower, unsigned mask){
  if (tower>15)return;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3;
  _towers[reg]&=(~(0x7<<bitshift))&0xfff;
  _towers[reg]|=(mask&0x7)<<bitshift;
}
const unsigned TrgDisabledChannels::towerMask(unsigned tower) const{
  if (tower>15)return 0;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3;
  return (_towers[reg]>>bitshift)&0x7;
}

void TrgDisabledChannels::setCalLowEnable(unsigned tower, bool on){
  if (tower>15)return;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3;
  if (on)_towers[reg]|=1<<bitshift;
  else _towers[reg]&=(~(1<<bitshift))&0xfff;
}
const bool TrgDisabledChannels::calLowEnabled(unsigned tower) const{
  if (tower>15)return 0;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3;
  return (_towers[reg]>>bitshift)&0x1;
}


void TrgDisabledChannels::setCalHighEnable(unsigned tower, bool on){
  if (tower>15)return;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3+1;
  if (on)_towers[reg]|=1<<bitshift;
  else _towers[reg]&=(~(1<<bitshift))&0xfff;
}
const bool TrgDisabledChannels::calHighEnabled(unsigned tower) const{
  if (tower>15)return 0;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3+1;
  return (_towers[reg]>>bitshift)&0x1;
}


void TrgDisabledChannels::setTkrEnable(unsigned tower, bool on){
  if (tower>15)return;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3+2;
  if (on)_towers[reg]|=1<<bitshift;
  else _towers[reg]&=(~(1<<bitshift))&0xfff;
}
const bool TrgDisabledChannels::tkrEnabled(unsigned tower) const{
  if (tower>15)return 0;
  unsigned reg=tower/4;
  unsigned bitshift=(tower%4)*3+2;
  return (_towers[reg]>>bitshift)&0x1;
}

void TrgDisabledChannels::setBusyFlag(unsigned tower, bool on){
  if (tower>15)return;
  if (on)_busy|=1<<tower;
  else _busy&=(~(1<<tower))&0xffff;
}
const bool TrgDisabledChannels::busyFlagEnabled(unsigned tower) const{
  if (tower>15)return 0;
  return (_busy>>tower)&0x1;
}

void TrgDisabledChannels::setTileEnableByIndex(unsigned tile, TrgConditions::pmt pmt,bool on){
  if (tile>108)return;
  unsigned reg=tile/9;
  unsigned bitshift=(tile%9);
  if(pmt ==TrgConditions::PMTB)bitshift+=9;
  if (on)_tiles[reg]|=1<<bitshift;
  else _tiles[reg]&=(~(1<<bitshift))&0x3ffff;
}
const bool TrgDisabledChannels::tileEnabledByIndex(unsigned tile, TrgConditions::pmt pmt) const{
  if (tile>108)return 0;
  unsigned reg=tile/9;
  unsigned bitshift=(tile%9);
  if(pmt ==TrgConditions::PMTB)bitshift+=9;
  return (_tiles[reg]>>bitshift)&0x1;
}

void TrgDisabledChannels::setTileEnableByName(unsigned tile, TrgConditions::pmt pmt,bool on){
  setTileEnableByIndex(AcdId::indexFromTile(tile),pmt,on);
}
const bool TrgDisabledChannels::tileEnabledByName(unsigned tile, TrgConditions::pmt pmt) const{
  return tileEnabledByIndex(AcdId::indexFromTile(tile),pmt);
}

// Attach this value to a TTree
void TrgDisabledChannels::makeBranch(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  
  std::string branchNameTower = branchName; branchNameTower += "towers";
  std::string leafNameTower = branchNameTower;
  leafNameTower += "[4]/i";
  tree.Branch(branchNameTower.c_str(),(void*)(_towers),leafNameTower.c_str());

  std::string branchNameTile = branchName; branchNameTile += "tiles";
  std::string leafNameTile = branchNameTile;
  leafNameTile += "[12]/i";
  tree.Branch(branchNameTile.c_str(),(void*)(_tiles),leafNameTile.c_str());

  std::string branchNameCno = branchName; branchNameCno += "acd_cno";
  std::string leafNameCno = branchNameCno;
  leafNameCno += "/i";
  tree.Branch(branchNameCno.c_str(),(void*)(&_cno),leafNameCno.c_str());

  std::string branchNameBusy = branchName; branchNameBusy += "tower_busy";
  std::string leafNameBusy = branchNameBusy;
  leafNameBusy += "/i";
  tree.Branch(branchNameBusy.c_str(),(void*)(&_busy),leafNameBusy.c_str());

  std::string branchNameExt = branchName; branchNameExt += "external";
  std::string leafNameExt = branchNameExt;
  leafNameExt += "/B";
  tree.Branch(branchNameExt.c_str(),(void*)(&_external),leafNameExt.c_str());

}

void TrgDisabledChannels::attach(TTree& tree, const std::string& prefix) const {
  std::string branchName = prefix; branchName += name();
  
  std::string branchNameTower = branchName; branchNameTower += "towers";
  tree.SetBranchAddress(branchNameTower.c_str(),(void*)(_towers));

  std::string branchNameTile = branchName; branchNameTile += "tiles";
  tree.SetBranchAddress(branchNameTile.c_str(),(void*)(_tiles));

  std::string branchNameCno = branchName; branchNameCno += "acd_cno";
  tree.SetBranchAddress(branchNameCno.c_str(),(void*)(&_cno));

  std::string branchNameBusy = branchName; branchNameBusy += "tower_busy";
  tree.SetBranchAddress(branchNameBusy.c_str(),(void*)(&_busy));

  std::string branchNameExt = branchName; branchNameExt += "external";
  tree.SetBranchAddress(branchNameExt.c_str(),(void*)(&_external));
}

std::ostream& operator <<(std::ostream& os, const TrgDisabledChannels* tc){
  os<<(*tc);
  return os;
} 
std::ostream& operator <<(std::ostream& os, const TrgDisabledChannels& tc){
  os<<"Disabled CAL low towers: ";
  for (int i=0;i<16;i++){
    if (!tc.calLowEnabled(i))os<<i<<" ";
  }
  os<<std::endl;
  os<<"Disabled CAL High towers: ";
  for (int i=0;i<16;i++){
    if (!tc.calHighEnabled(i))os<<i<<" ";
  }
  os<<std::endl;
  os<<"Disabled TKR towers: ";
  for (int i=0;i<16;i++){
    if (!tc.tkrEnabled(i))os<<i<<" ";
  }
  os<<std::endl;
  os<<"Disabled ACD tiles: ";
  for (int i=0;i<108;i++){
    unsigned tile=AcdId::tileFromIndex(i);
    if (!tc.tileEnabledByName(tile,TrgConditions::PMTA))os<<std::setfill('0')<<std::setw(3)<<tile<<"-"<<"A"<<" ";
    if (!tc.tileEnabledByName(tile,TrgConditions::PMTB))os<<std::setfill('0')<<std::setw(3)<<tile<<"-"<<"B"<<" ";
  }
  os<<std::endl;
  os<<"Disabled Busy signals: ";
  for (int i=0;i<16;i++){
    if (!tc.busyFlagEnabled(i))os<<i<<" ";
  }
  os<<std::endl;
  if (tc.externalEnabled())os<<"External enabled";
  else os<<"External disabled";
  os<<std::endl;
  return os;
}
      
      


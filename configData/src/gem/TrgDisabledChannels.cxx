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

TrgDisabledChannels::TrgDisabledChannels(){
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
  _busy=val&0xff;
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
      
      


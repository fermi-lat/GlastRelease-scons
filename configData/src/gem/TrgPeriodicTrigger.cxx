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

TrgPeriodicTrigger::TrgPeriodicTrigger(){
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

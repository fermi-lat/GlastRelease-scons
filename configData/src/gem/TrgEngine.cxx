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


TrgEngine::TrgEngine(){
  _engine=0;
}
TrgEngine::TrgEngine(unsigned long en):_engine(en){ }

void TrgEngine::setEngine(unsigned long en){
  _engine=en;
}
unsigned long TrgEngine::engine(){
  return _engine;
}
void TrgEngine::setCalstrobe(bool on){
  if(on)_engine|=1<<16;
  else _engine&=~(1<<16);
}
void TrgEngine::setTack(bool on){
  if(on)_engine|=1<<19;
  else _engine&=~(1<<19);
}
void TrgEngine::setZeroSuppression(bool on){
  if(on)_engine|=1<<21;
  else _engine&=~(1<<21);
}
void TrgEngine::setFourRangeReadout(bool on){
  if(on)_engine|=1<<20;
  else _engine&=~(1<<20);

}
void TrgEngine::setInhibit(bool on){
  if(on)_engine|=1<<15;
  else _engine&=~(1<<15);
}
void TrgEngine::setDestination(int dest){
  _engine&=~(0x3f<<25);
  _engine|=((dest&0x3f)<<25);
}
void TrgEngine::setMarker(int marker){
  _engine&=~(0x7<<22);
  _engine|=((marker&0x7)<<22);
}
void TrgEngine::setPrescale(int num){
  _engine&=~0xff;
  _engine|=num&0xff;
}
bool TrgEngine::calstrobe() const{
  return _engine&(1<<16)? 1: 0;
}
bool TrgEngine::tack() const{
  return _engine&(1<<19)? 1: 0;
}
bool TrgEngine::inhibited() const{
  return _engine&(1<<15)? 1: 0;
}
bool TrgEngine::zeroSuppression() const{
  return _engine&(1<<21)? 1: 0;
}
bool TrgEngine::fourRangeReadout() const{
  return _engine&(1<<20)? 1: 0;
}
int TrgEngine::destination() const{
  return (_engine>>25)&0x3f;
}
int TrgEngine::marker() const{
  return (_engine>>22)&0x7;
}
int TrgEngine::prescale() const{
  return _engine&0xff;
}
void TrgEngine::clear(){
  _engine=0;
} 

std::ostream& operator <<(std::ostream& os, const TrgEngine& tc){
  if (tc.calstrobe())os<<"Calstrobe ";
  if (tc.tack())os<<"Tack ";
  if (tc.inhibited())os<<"Inhibit ";
  os<<"; ";
  if (tc.prescale()!=0)os<<"Prescale="<<tc.prescale()<<" ";
  if (tc.marker()!=0)os<<"Marker="<<tc.marker()<<" ";
  if (tc.destination()!=0)os<<"Destination="<<tc.destination()<<" ";
  if (tc.zeroSuppression())os<<"0-suppression ";
  if (tc.fourRangeReadout())os<<"4-range-readout ";
  return os;
}
    


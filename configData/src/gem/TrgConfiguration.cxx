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

TrgConfiguration::TrgConfiguration():_configuration(0){}

TrgConfiguration::TrgConfiguration(unsigned long c):_configuration(c){}

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

std::ostream& operator <<(std::ostream& os,  const TrgConfiguration& tc){
  if(tc.acdUsedAsTrigger())os<<"ACD used as trigger"<<std::endl;
  else os<<"ACD used as veto"<<std::endl;
  os<<"Trigger word "<<std::hex<<tc.configuration()<<std::dec<<std::endl;
  return os;
}

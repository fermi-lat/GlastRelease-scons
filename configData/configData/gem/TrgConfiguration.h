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

#ifndef TRGCONFIGURATION_HH
#define TRGCONFIGURATION_HH

#include <iostream>

class TrgConfiguration{
 public:
  TrgConfiguration();
  TrgConfiguration(unsigned long i);
  void setConfiguration(unsigned long i);
  unsigned long configuration() const{return _configuration;}
  void useAcdAsTrigger(bool on);
  bool acdUsedAsTrigger() const;
 private:
  unsigned long _configuration;
};
#endif
std::ostream& operator <<(std::ostream& os, const TrgConfiguration& tc);

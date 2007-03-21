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
#include "configData/base/ConfigTuple.h"

class TrgConfiguration : public ConfigBranch {
public:
  TrgConfiguration();
  TrgConfiguration(unsigned long i);
  virtual ~TrgConfiguration(){}
  void setConfiguration(unsigned long i);
  unsigned long configuration() const{return _configuration;}
  void useAcdAsTrigger(bool on);
  bool acdUsedAsTrigger() const;
  // Reset the cached and output values
  virtual void reset(){_configuration=0;}
  // Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;
private:
  unsigned long _configuration;
};
#endif
std::ostream& operator <<(std::ostream& os, const TrgConfiguration& tc);

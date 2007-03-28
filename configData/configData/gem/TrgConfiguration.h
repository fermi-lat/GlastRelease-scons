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
/** @class TrgConfiguration
    @author M. Kocian

    Class to hold GEM configuration register information.
*/


#ifndef TRGCONFIGURATION_HH
#define TRGCONFIGURATION_HH

#include <iostream>
#include "configData/base/ConfigTuple.h"

class TrgConfiguration : public ConfigBranch {
public:
  /// Default constructor
  TrgConfiguration();
  /// Constructor with a register value as initializer
  TrgConfiguration(unsigned long i);
  /// Default destructor
  virtual ~TrgConfiguration(){}
  /// Set configuration to register value i
  void setConfiguration(unsigned long i);
  /// Retrieve configuration register value
  unsigned long configuration() const{return _configuration;}
  /// Set bit so ACD is used as trigger 
  void useAcdAsTrigger(bool on);
  /// Is ACD used as trigger?
  bool acdUsedAsTrigger() const;
  /// Reset the cached and output values
  virtual void reset(){_configuration=0;}
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;
private:
  unsigned long _configuration;
};
#endif
/// printout function
std::ostream& operator <<(std::ostream& os, const TrgConfiguration& tc);

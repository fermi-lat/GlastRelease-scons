//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Trigger window parameters in GEM configuration
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
/** @class TrgWindowParams
    @author M. Kocian

    Class to hold GEM window parameters.
*/
#ifndef TRGWINDOWPARAMS_HH
#define TRGWINDOWPARAMS_HH

#include "TrgConditions.h"
#include <iostream>
#include "configData/base/ConfigTuple.h"

class TrgWindowParams : public ConfigBranch {
  
public:
  /// Default constructor
  TrgWindowParams();
  /// Default destructor
  virtual ~TrgWindowParams(){}
  /// Retrieve window mask
  int windowMask() const;
  /// Retrieve window width
  int windowWidth() const;
  /// Query if the mask includes a particular condition
  bool maskIncludes(TrgConditions::conditions);
  /// Reset the cached and output values
  virtual void reset(){_windowmask=0;_winwidth=0;}
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;

private:
  friend class TrgConfigParser;

  /// set window mask to a value
  void setWindowMask(int);
  /// set window width to a value
  void setWindowWidth(int);
  /// add a condition to the window open mask
  void addToOpenMask(TrgConditions::conditions);

  char _windowmask;
  int _winwidth;
};
/// printout function  
std::ostream& operator <<(std::ostream &os, const TrgWindowParams& tl);

#endif

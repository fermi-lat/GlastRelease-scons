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

#ifndef TRGWINDOWPARAMS_HH
#define TRGWINDOWPARAMS_HH

#include "TrgConditions.h"
#include <iostream>
#include "configData/base/ConfigTuple.h"

class TrgWindowParams : public ConfigBranch {
  
public:
  TrgWindowParams();
  virtual ~TrgWindowParams(){}
  int windowMask() const;
  int windowWidth() const;
  bool maskIncludes(TrgConditions::conditions);
  // Reset the cached and output values
  virtual void reset(){_windowmask=0;_winwidth=0;}
  // Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;

private:
  friend class TrgConfigParser;

  void setWindowMask(int);
  void setWindowWidth(int);
  void addToOpenMask(TrgConditions::conditions);

  char _windowmask;
  int _winwidth;
};
  
std::ostream& operator <<(std::ostream &os, const TrgWindowParams& tl);

#endif

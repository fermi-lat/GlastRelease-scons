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

class TrgWindowParams{
  
public:
  TrgWindowParams();
  int windowMask() const;
  int windowWidth() const;
  bool maskIncludes(TrgConditions::conditions);
  
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

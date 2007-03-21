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
#ifndef TRGENGINE_HH
#define TRGENGINE_HH

#include <iostream>
#include "configData/base/ConfigTuple.h"

class TrgEngine : public ConfigBranch {

public:
  TrgEngine();
  void setEngine(unsigned idx, unsigned long engine);
  unsigned long engine(unsigned idx) const;
  void setCalstrobe(unsigned idx, bool);
  bool calstrobe(unsigned idx) const;
  void setTack(unsigned idx, bool);
  bool tack(unsigned idx) const;
  void setZeroSuppression(unsigned idx, bool);
  bool zeroSuppression(unsigned idx) const;
  void setFourRangeReadout(unsigned idx, bool);
  bool fourRangeReadout(unsigned idx) const;
  void setInhibit(unsigned idx, bool);
  bool inhibited(unsigned idx) const;
  void setDestination(unsigned idx, int);
  int destination(unsigned idx) const;
  void setMarker(unsigned idx, int);
  int marker(unsigned idx) const;
  int prescale(unsigned idx) const;
  void setPrescale(unsigned idx, int);
  void clear();

  // Reset the cached and output values
  virtual void reset() { clear() ;}
  // Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;    

private:

  unsigned long _engine[16];
};

std::ostream& operator <<(std::ostream& os, const TrgEngine& tc);

#endif

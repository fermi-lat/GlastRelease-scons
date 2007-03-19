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

class TrgEngine{

public:
  TrgEngine();
  TrgEngine(unsigned long engine);
  void setEngine(unsigned long engine);
  unsigned long engine();
  void setCalstrobe(bool);
  bool calstrobe() const;
  void setTack(bool);
  bool tack() const;
  void setZeroSuppression(bool);
  bool zeroSuppression() const;
  void setFourRangeReadout(bool);
  bool fourRangeReadout() const;
  void setInhibit(bool);
  bool inhibited() const;
  void setDestination(int);
  int destination() const;
  void setMarker(int);
  int marker() const;
  int prescale() const;
  void setPrescale(int);
  void clear();
  

private:

  unsigned long _engine;
};

std::ostream& operator <<(std::ostream& os, const TrgEngine& tc);

#endif

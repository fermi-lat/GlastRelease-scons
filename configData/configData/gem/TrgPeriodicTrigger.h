//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Periodic Trigger Parameters in GEM configuration
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

#ifndef TRGPERIODICTRIGGER_HH
#define TRGPERIODICTRIGGER_HH
#include <iostream>

class TrgPeriodicTrigger{
public:
  TrgPeriodicTrigger();
  int prescale() const;
  bool onePPS() const;
  bool freeRunning() const;
  int limit() const;

private:

  void setPrescale(int);
  void setOnePPS(bool);
  void setFreeRunning(bool);
  void setLimit(int);
  
  int _prescale;
  int _limit;
  bool _free;
  bool _onepps;

  friend class TrgConfigParser;
};

std::ostream& operator <<(std::ostream &os, const TrgPeriodicTrigger& tl);

#endif

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
/** @class TrgPeriodicTrigger
    @author M. Kocian

    Class to hold GEM periodic trigger information.
*/

#ifndef TRGPERIODICTRIGGER_HH
#define TRGPERIODICTRIGGER_HH
#include <iostream>
#include "configData/base/ConfigTuple.h"

class TrgPeriodicTrigger : public ConfigBranch {
public:
  /// Default constructor
  TrgPeriodicTrigger();
  /// return the prescale (which determines the frequency)
  int prescale() const;
  /// Is the basic time unit 1-PPS?
  bool onePPS() const;
  /// Is the trigger free running? This is actually controlled by FSW so 
  /// it may not be a good idea to use this function
  bool freeRunning() const;
  /// Number of triggers for limited mode
  int limit() const;
  /// Reset the cached and output values
  virtual void reset();
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;

private:

  /// Set prescale (which determines the frequency)
  void setPrescale(int);
  /// Set basic time unit to 1-PPS (true) or internal clock=50 ns (false)
  void setOnePPS(bool);
  /// Set the periodic trigger to be free running or not.
  /// This is really controlled by FSW
  void setFreeRunning(bool);
  /// Set the number of triggers in limited mode.
  void setLimit(int);
  
  int _prescale;
  int _limit;
  bool _free;
  bool _onepps;

  friend class TrgConfigParser;
};

/// printout function
std::ostream& operator <<(std::ostream &os, const TrgPeriodicTrigger& tl);

#endif

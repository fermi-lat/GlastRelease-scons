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
/** @class TrgEngine
    @author M. Kocian

    Class to hold information about the 16 trigger engines of the GEM
*/

class TrgEngine : public ConfigBranch {

public:
  /// Default constructor
  TrgEngine();
  /// Set value of the register for engine idx to value engine
  void setEngine(unsigned idx, unsigned long engine);
  /// Retrieve value for the register of engine idx
  unsigned long engine(unsigned idx) const;
  /// Set calstrobe bit of engine idx
  void setCalstrobe(unsigned idx, bool);
  /// Retrieve calstrobe bit of engine idx
  bool calstrobe(unsigned idx) const;
  /// Set TACK bit of engine idx
  void setTack(unsigned idx, bool);
  /// Retrieve TACK bit of engine idx
  bool tack(unsigned idx) const;
  /// Set 0 suppression bit of engine idx
  void setZeroSuppression(unsigned idx, bool);
  /// Retrieve 0 suppression bit of engine idx
  bool zeroSuppression(unsigned idx) const;
  /// Set 4-range readout bit of engine idx
  void setFourRangeReadout(unsigned idx, bool);
  /// Retrieve 4-range readout bit of engine idx
  bool fourRangeReadout(unsigned idx) const;
  /// Set inhibit bit of engine idx
  void setInhibit(unsigned idx, bool);
  /// Retrieve inhibit bit of engine idx
  bool inhibited(unsigned idx) const;
  /// Set destination of engine idx
  void setDestination(unsigned idx, int);
  /// Retrieve destination of engine idx
  int destination(unsigned idx) const;
  /// Set marker of engine idx
  void setMarker(unsigned idx, int);
  /// Retrieve marker of engine idx
  int marker(unsigned idx) const;
  /// Retrieve prescale of engine idx
  int prescale(unsigned idx) const;
  /// Set prescale of engine idx
  void setPrescale(unsigned idx, int);
  /// Clear data
  void clear();

  /// Reset the cached and output values
  virtual void reset() { clear() ;}
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;    

private:

  unsigned long _engine[16];
};

/// printout function
std::ostream& operator <<(std::ostream& os, const TrgEngine& tc);

#endif

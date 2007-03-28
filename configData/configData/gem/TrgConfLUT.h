//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      GEM LUT class
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
#ifndef TRGCONFLUT_HH
#define TRGCONFLUT_HH
#include <string.h>
#include <iostream>
#include <vector>
#include "TrgConditions.h"

#include "configData/base/ConfigTuple.h"
/** @class TrgConfLUT
    @author M. Kocian

    Class to hold the GEM lookup table that links trigger conditions
    to trigger engines.
*/

// number of registers in LUT
#define NREG 32
// number of conditions per register
#define CONDSREG 8
// number of entries in LUT
#define NENTRIES NREG*CONDSREG
// number of bits per condition
#define NBITS 4
// mask one condition (NBITS bits high)
#define CONDMASK 0xf
//---------------

class TrgConfLUT: public ConfigBranch {

public:
  /// Default constructor
  TrgConfLUT();
  /// Returns the engine number for condition i
  int engineNumber(int i) const;
  /// Returns a vector of conditions for engine i
  std::vector<unsigned long> conditionsInEngine(int i) const;
  /// Returns the content of the ith CONDITIONS register in the GEM
  unsigned long reg(int i) const {if (i<NREG && i>=0)return _lut[i];
                                       else return 0;} 

  /// Reset the cached and output values
  virtual void reset();
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;    

private:
  friend class TrgConfigParser;
  /// Point condition i to engine entry
  void setEngineNumber(int i, int entry);
  /// Set value corresponding to the ith GEM CONDITIONS register to value r
  void setRegister(int i, unsigned long r){if (i<NREG && i>=0)_lut[i]=r;}

  unsigned long _lut[NREG];
};
std::ostream& operator <<(std::ostream& os, const TrgConfLUT& tc);
#endif

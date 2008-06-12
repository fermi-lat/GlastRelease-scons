//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      A GEM ROI register class 
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
/** @class FswEfcSampler
    @author M. Kocian

    Class to hold GEM ROI information.
*/

#ifndef FSWEFCSAMPLER_HH
#define FSWEFCSAMPLER_HH

#include <iostream>
#include <vector>

#include "configData/base/ConfigTuple.h"
#include "enums/Lsf.h"

class FswEfcSampler: public ConfigBranch {

public:

  static FswEfcSampler* makeFromXmlFile(const char* fileName);

  static unsigned mapKey(unsigned lpaMode, unsigned handlerId);

public:
  /// Default constructor
  FswEfcSampler();
  /// Clear data contents
  void clear();

  /// 
  inline unsigned prescaler(unsigned line) const { return m_prescalers[line]; }
  ///
  inline unsigned inputPrescaler() const { return m_inputPrescaler; }
  ///
  inline unsigned outputPrescaler() const { return m_outputPrescaler; }
  ///
  inline unsigned enabled() const { return m_enabled; }
  ///

  void set(unsigned prescalers[32], unsigned input, unsigned output, unsigned enabled);

  // Return the prescale factor for a given event
  unsigned prescaleFactor(enums::Lsf::RsdState rsdState, enums::Lsf::LeakedPrescaler leakedPrescaler) const;


  /// Reset the cached and output values
  virtual void reset() { clear() ;}
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;  
private:
  unsigned m_prescalers[32];
  unsigned m_inputPrescaler;
  unsigned m_outputPrescaler;
  unsigned m_enabled;
}; 
#endif
/// Printout Function
std::ostream& operator <<(std::ostream& os, const FswEfcSampler& tc);

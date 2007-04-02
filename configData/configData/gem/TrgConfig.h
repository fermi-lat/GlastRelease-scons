//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      A top level GEM configuration class 
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
/** @class TrgConfig
    @author M. Kocian

    Top level class to hold GEM information.
*/

#ifndef TRGCONFIG_HH
#define TRGCONFIG_HH

#include "TrgConfLUT.h"
#include "TrgEngine.h"
#include "TrgDisabledChannels.h"
#include "TrgWindowParams.h"
#include "TrgPeriodicTrigger.h"
#include "TrgConfiguration.h"
#include "TrgRoi.h"

class TrgConfig{

public:
  /// Default constructor
  TrgConfig();
  /// Return the window parameter object
  TrgWindowParams* windowParams(){return &_twp;}
  /// Return the periodic trigger object
  TrgPeriodicTrigger* periodicTrigger(){return &_tpt;}
  /// Return the lookup table object
  TrgConfLUT* lut(){return &_lut;}
  /// Return the trigger engine object
  TrgEngine* trgEngine() {return &_tev;}
  /// Return the object that holds the information about input disables
  TrgDisabledChannels* disabledChannels(){return &_tdv;}
  /// Return the information contained in the configuration register
  TrgConfiguration* configuration(){return &_configuration;}
  /// Return the object that holds the ROI configuration
  TrgRoi* roi(){return &_roi;}

  /// Reset the cached and output values
  void reset();
  /// Attach this value to a TTree
  void makeBranches(TTree& tree, const std::string& prefix) const;
  void attach(TTree& tree, const std::string& prefix) const;
  
private:
  TrgWindowParams _twp;
  TrgPeriodicTrigger _tpt;
  TrgConfLUT _lut;
  TrgEngine _tev;
  TrgDisabledChannels _tdv;
  TrgConfiguration _configuration;
  TrgRoi _roi;
};

/// printout function  
std::ostream& operator <<(std::ostream &os, const TrgConfig& tl);




#endif

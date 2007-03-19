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
  TrgConfig();
  TrgWindowParams* windowParams(){return &_twp;}
  TrgPeriodicTrigger* periodicTrigger(){return &_tpt;}
  TrgConfLUT* lut(){return &_lut;}
  TrgEngine* trgEngine(int i);
  TrgDisabledChannels* disabledChannels(){return &_tdv;}
  TrgConfiguration* configuration(){return &_configuration;}
  TrgRoi* roi(){return &_roi;}
  
private:
  TrgWindowParams _twp;
  TrgPeriodicTrigger _tpt;
  TrgConfLUT _lut;
  std::vector<TrgEngine> _tev;
  TrgDisabledChannels _tdv;
  TrgConfiguration _configuration;
  TrgRoi _roi;
};





#endif

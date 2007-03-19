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
#include "configData/gem/TrgConfig.h"


TrgConfig::TrgConfig():_tev(16){}

TrgEngine* TrgConfig::trgEngine(int i){
  if (i<16)return &_tev[i];
  else return 0;
}

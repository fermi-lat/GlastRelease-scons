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

#ifndef TRGROI_HH
#define TRGROI_HH

#include <iostream>
#include <vector>

#include "configData/base/ConfigTuple.h"

class TrgRoi: public ConfigBranch {
public:
  TrgRoi();
  void clear();
  void setRoiRegister(int i, unsigned long val);
  unsigned long roiRegister(int i) const;
  std::vector<unsigned long> roiFromIndex(int index) const;
  std::vector<unsigned long> roiFromName(int name) const;
  std::vector<unsigned long> tilesByIndexInRoi(int i) const;
  std::vector<unsigned long> tilesByNameInRoi(int i) const;
  bool isTileByIndexInRoi(unsigned tile, int roi) const;
  bool isTileByNameInRoi(unsigned tile, int roi) const;
  // Reset the cached and output values
  virtual void reset() { clear() ;}
  // Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;  
private:
  unsigned long _roi[54];
};
#endif
std::ostream& operator <<(std::ostream& os, const TrgRoi& tc);

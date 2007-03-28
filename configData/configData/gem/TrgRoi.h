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
/** @class TrgRoi
    @author M. Kocian

    Class to hold GEM ROI information.
*/

#ifndef TRGROI_HH
#define TRGROI_HH

#include <iostream>
#include <vector>

#include "configData/base/ConfigTuple.h"

class TrgRoi: public ConfigBranch {
public:
  /// Default constructor
  TrgRoi();
  /// Clear data contents
  void clear();
  /// Set ROI register number i to value val
  void setRoiRegister(int i, unsigned long val);
  /// Retrieve value of ROI register number i
  unsigned long roiRegister(int i) const;
  /// Retrieve vector of ROIs that contain tile index 
  /// where index refers to tile number as defined in TD-01545 Appendix A
  std::vector<unsigned long> roiFromIndex(int index) const;
  /// Retrieve vector of ROIs that contain tile name
  /// where name is the tile name (doh!)
  std::vector<unsigned long> roiFromName(int name) const;
  /// Retrieve vector of tile indeces for ROI i
  /// where index refers to tile number as defined in TD-01545 Appendix A
  std::vector<unsigned long> tilesByIndexInRoi(int i) const;
  /// Retrieve vector of tile names for ROI i
  std::vector<unsigned long> tilesByNameInRoi(int i) const;
  /// Is tile with index tile in roi roi?
  /// Index refers to tile number as defined in TD-01545 Appendix A
  bool isTileByIndexInRoi(unsigned tile, int roi) const;
  /// Is tile with name tile in roi roi?
  bool isTileByNameInRoi(unsigned tile, int roi) const;
  /// Reset the cached and output values
  virtual void reset() { clear() ;}
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;  
private:
  unsigned long _roi[54];
}; 
#endif
/// Printout Function
std::ostream& operator <<(std::ostream& os, const TrgRoi& tc);

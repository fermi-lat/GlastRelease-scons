//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Disabled channels in GEM configuration
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

#ifndef TRGDISABLEDCHANNELS_HH
#define TRGDISABLEDCHANNELS_HH

#include <vector>
#include <iostream>
#include "TrgConditions.h"
#include "configData/base/ConfigTuple.h"

/** @class TrgDisabledChannels
    @author M. Kocian

    Class to hold information about disabled inputs in the GEM
*/

class TrgDisabledChannels : public ConfigBranch {

public:
  /// Default constructor
  TrgDisabledChannels();
  /// clear data
  void clear();
  /// Set ith TOWERS GEM register to value val
  void setTowerRegister(unsigned i,unsigned val);
  /// Retrieve ith TOWERS GEM register
  unsigned towerRegister(unsigned i);
  /// Set ith TILES GEM register to value val
  void setTileRegister(unsigned i, unsigned val);
  /// Retrieve ith TILES GEM register
  unsigned tileRegister(unsigned i);
  /// Set CNO GEM register to value val
  void setCnoRegister(unsigned val);
  /// Retrieve CNO GEM register
  unsigned cnoRegister(){return _cno;}
  /// Set TOWER_BUSY GEM register to value val
  void setBusyRegister(unsigned val);
  /// Retrieve TOWER_BUSY GEM register 
  unsigned busyRegister(){return _busy;}
  /// Set EXTERNAL GEM register to value val
  void setExternal(unsigned val);
  /// Retrieve EXTERNAL GEM register 
  const bool externalEnabled() const{return _external;}
  /// Set enable mask for tower tower to mask 
  void setTowerMask(unsigned tower, unsigned mask);
  /// Retrieve enable mask for tower tower
  const unsigned towerMask(unsigned tower) const;
  /// Toggle CALLO enable for tower tower
  void setCalLowEnable(unsigned tower, bool on);
  /// Toggle CALHI enable for tower tower
  void setCalHighEnable(unsigned tower, bool on);
  /// Toggle TKR enable for tower tower
  void setTkrEnable(unsigned tower, bool on);
  /// Retrieve CALLO enable bit for tower tower
  const bool calLowEnabled(unsigned tower) const;
  /// Retrieve CALHI enable bit for tower tower
  const bool calHighEnabled(unsigned tower) const;
  /// Retrieve TKR enable bit for tower tower
  const bool tkrEnabled(unsigned tower) const;
  /// Toggle Busy flag for tower tower
  void setBusyFlag(unsigned tower, bool on);
  /// Retrieve Busy flag for tower tower
  const bool busyFlagEnabled(unsigned tower) const;
  /// Toggle tile enable for the tile with index tile and PMT pmt
  /// (This is the index which is called Tile Number defined in Appendix A of TD-01545,
  /// not the index as used with the GEM event contribution)
  void setTileEnableByIndex(unsigned tile, TrgConditions::pmt pmt, bool on);
  /// Retrieve tile enable for the tile with index tile and PMT pmt
  /// (This is the index which is called Tile Number defined in Appendix A of TD-01545,
  /// not the index as used with the GEM event contribution)
  const bool tileEnabledByIndex(unsigned tile,TrgConditions::pmt pmt) const;
  /// Toggle tile enable for the tile with name tile and PMT pmt
  void setTileEnableByName(unsigned tile, TrgConditions::pmt pmt, bool on);
  /// Retrieve tile enable for the tile with name tile and PMT pmt
  const bool tileEnabledByName(unsigned tile,TrgConditions::pmt pmt) const;
  
  /// Reset the cached and output values
  virtual void reset() { clear() ;}
  /// Attach this value to a TTree
  virtual void makeBranch(TTree& tree, const std::string& prefix) const;
  virtual void attach(TTree& tree, const std::string& prefix) const;

private:
  
  unsigned _towers[4];
  unsigned _tiles[12];
  unsigned _cno;
  unsigned _busy;
  bool _external;

  friend class TrgConfigParser;
};
/// printout function
std::ostream& operator <<(std::ostream& os, const TrgDisabledChannels& tc);
/// printout function
std::ostream& operator <<(std::ostream& os, const TrgDisabledChannels* tc);


#endif

/**
   @file TkrSplitCalib.cxx

   Implementation file for calib. TDS class for Tracker splits

   @author J. Bogart
*/
#include "CalibData/Tkr/TkrSplitsCalib.h"

namespace CalibData {


  void TkrSplit::update(RangeBase* other) {

    // had better work!
    TkrSplit* otherSplit = dynamic_cast<TkrSplit* > (other);

    m_low = otherSplit->m_low;
    m_high = otherSplit->m_high;
  }

  const CLID&  TkrSplitsCalib::classID() {return CLID_Calib_TKR_Splits;}

  RangeBase* TkrSplitsCalib::getChannel(const idents::TkrId& id, 
                                        unsigned feChip) { 
    // Following sets index 
    TkrBase::getChannel(id, feChip);
    if (!m_ixValid) return 0;

    return &m_splits[m_ix];
}
    
  RangeBase* TkrSplitsCalib::getChannel(unsigned towerRow, unsigned towerCol,
                                        unsigned tray, bool top, 
                                        unsigned feChip) {
    idents::TkrId id(towerCol, towerRow, tray, top);
    return getChannel(id, feChip);
  }
    
    
  bool TkrSplitsCalib::putChannel(RangeBase* data, const idents::TkrId& id, 
                                  unsigned feChip) {
    bool ok = TkrBase::putChannel(data, id, feChip);
    if (!ok) return false;

    m_splits[m_ix].update(data);
    return true;
  }


  bool TkrSplitsCalib::putChannel(RangeBase* data, unsigned towerRow, 
                                  unsigned towerCol, unsigned tray, 
                                  bool top, unsigned feChip) {
    idents::TkrId id(towerCol, towerRow, tray, top);
    return putChannel(data, id, feChip);
  }


  StatusCode TkrSplitsCalib::update(CalibBase& other, MsgStream* log){
    StatusCode ret = TkrBase::update(other, log);
    if (!ret.isSuccess()) return ret;

    TkrSplitsCalib& other1 = dynamic_cast<TkrSplitsCalib& >(other);

    unsigned n = m_finder->getSize();
    for (unsigned i = 0; i < n; i++) {
      m_splits[i].update(&(other1.m_splits[i]) );
    }
    return StatusCode::SUCCESS;
  }

}



// $Header$

#include "CalibData/Cal/CalCalibLightAsym.h"
#include "CalibData/Cal/CalFinder.h"
#include "CalibData/CalibModel.h"

namespace {
  unsigned diode(unsigned range) {
    if ((range == idents::CalXtalId::LEX8) || 
        (range == idents::CalXtalId::LEX1) ) return 0;
    return 1;
  }
}
namespace CalibData {

  CalCalibLightAsym::CalCalibLightAsym(unsigned nTowerRow, unsigned nTowerCol, 
                             unsigned nLayer, unsigned nXtal, 
                             unsigned nFace, unsigned nRange) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, nFace, nRange),
    m_nRange(0) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();
    m_nRange = m_finder->getNRange();
    
    LightAsym* pLightAsyms = new LightAsym[size];
    for (ix = 0; ix < size; ix++) {
      m_ranges[ix] = pLightAsyms; 
      ++pLightAsyms;
    }
  }

  CalCalibLightAsym::~CalCalibLightAsym() {
    LightAsym* pLightAsyms = dynamic_cast<LightAsym* >(m_ranges[0]);
    if (pLightAsyms) delete [] pLightAsyms;
  }

  const CLID& CalCalibLightAsym::classID()   {return CLID_Calib_CAL_LightAsym;}

  RangeBase* CalCalibLightAsym::getRange(idents::CalXtalId id, unsigned range, 
                                         unsigned face)
  {
    if (m_nRange == 2) return CalCalibBase::getRange(id, diode(range), face);
    else return CalCalibBase::getRange(id, range, face);
  }

  bool CalCalibLightAsym::putRange(idents::CalXtalId id, unsigned range, 
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<LightAsym* >(data)) return false;

    // For Light asym, don't have separate data for all 4 range. Have
    // only two sets: one for LEX1, LEX8 and another for HEX1, HEX8
    // So map LEX1, LEX8 to 0; HEX1, HEX8 to 1

    // Otherwise go ahead and let base class handle it
    if (m_nRange == 2) {
      return CalCalibBase::putRange(id, diode(range), face, data);
    }
    else {
      return CalCalibBase::putRange(id, range, face, data);
    }
  }
  bool CalCalibLightAsym::putRange(unsigned towerRow, unsigned towerCol, 
                              unsigned layer, unsigned xtal, unsigned range,
                              unsigned face, RangeBase* data) {
    if (!dynamic_cast<LightAsym* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    if (m_nRange == 2) {
      return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                    diode(range), face, data);
    } else {
      return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                    range, face, data);
    }
    
  }

}

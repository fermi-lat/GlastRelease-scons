// $Header$

#include "CalibData/Cal/CalCalibMuSlope.h"
#include "CalibData/Cal/CalFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  CalCalibMuSlope::CalCalibMuSlope(unsigned nTowerRow, unsigned nTowerCol, 
                             unsigned nLayer, unsigned nXtal, 
                             unsigned nFace, unsigned nRange) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, nFace, nRange) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    MuSlope* pMuSlopes = new MuSlope[size];
    for (ix = 0; ix < size; ix++) {
      m_ranges[ix] = pMuSlopes; 
      ++pMuSlopes;
    }
  }

  CalCalibMuSlope::~CalCalibMuSlope() {
    MuSlope* pMuSlopes = dynamic_cast<MuSlope* >(m_ranges[0]);
    if (pMuSlopes) delete [] pMuSlopes;
  }

  const CLID& CalCalibMuSlope::classID()   {return CLID_Calib_CAL_MuSlope;}

  bool CalCalibMuSlope::putRange(idents::CalXtalId id, unsigned range, 
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<MuSlope* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }
  bool CalCalibMuSlope::putRange(unsigned towerRow, unsigned towerCol, 
                              unsigned layer, unsigned xtal, unsigned range,
                              unsigned face, RangeBase* data) {
    if (!dynamic_cast<MuSlope* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);
  }

}

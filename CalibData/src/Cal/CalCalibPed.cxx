// $Header$

#include "CalibData/Cal/CalCalibPed.h"
#include "CalibData/Cal/CalFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  CalCalibPed::CalCalibPed(unsigned nTowerRow, unsigned nTowerCol, 
                           unsigned nLayer, unsigned nXtal, 
                           unsigned nFace, unsigned nRange) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, nFace, nRange) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    for (ix = 0; ix < size; ix++) {
      (*m_pR)[ix] = new Ped(0.0, 0.0);
    }
  }


  const CLID& CalCalibPed::classID()   {return CLID_Calib_CAL_Ped;}

  bool CalCalibPed::putRange(idents::CalXtalId id, unsigned range, 
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<Ped* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

  bool CalCalibPed::putRange(unsigned towerRow, unsigned towerCol, 
                             unsigned layer, unsigned xtal, unsigned range,
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<Ped* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);
  }

}

// $Header$

#include "CalibData/Acd/AcdCalibPed.h"
#include "CalibData/Acd/AcdFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AcdCalibPed::AcdCalibPed(unsigned nFace, unsigned nRow, unsigned nCol, 
                             unsigned nPmt, unsigned nRange) :
    AcdCalibBase(nFace, nRow, nCol, nPmt, nRange) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AcdPed* pPeds = new AcdPed[size];
    for (ix = 0; ix < size; ix++) {
      m_ranges[ix] = pPeds; 
      ++pPeds;
    }
  }

  AcdCalibPed::~AcdCalibPed() {
    AcdPed* pPeds = dynamic_cast<AcdPed* >(m_ranges[0]);
    if (pPeds) delete [] pPeds;
  }

  const CLID& AcdCalibPed::classID()   {return CLID_Calib_ACD_Ped;}

  bool AcdCalibPed::putRange(idents::AcdId id, unsigned pmt, unsigned range, 
                              RangeBase* data) {
    if (!dynamic_cast<AcdPed* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putRange(id, pmt, range, data);
  }
  bool AcdCalibPed::putRange(unsigned face, unsigned row, unsigned col, 
                              unsigned pmt, unsigned range,
                              RangeBase* data) {
    if (!dynamic_cast<AcdPed* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putRange(face, row, col, pmt, range, data);
  }

}

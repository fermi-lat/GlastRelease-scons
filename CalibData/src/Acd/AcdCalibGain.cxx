// $Header$

#include "CalibData/Acd/AcdCalibGain.h"
#include "CalibData/Acd/AcdFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AcdCalibGain::AcdCalibGain(unsigned nFace, unsigned nRow, unsigned nCol, 
                             unsigned nPmt, unsigned nRange) :
    AcdCalibBase(nFace, nRow, nCol, nPmt, nRange) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AcdGain* pGains = new AcdGain[size];
    for (ix = 0; ix < size; ix++) {
      m_ranges[ix] = pGains; 
      ++pGains;
    }
  }

  AcdCalibGain::~AcdCalibGain() {
    AcdGain* pGains = dynamic_cast<AcdGain* >(m_ranges[0]);
    if (pGains) delete [] pGains;
  }

  const CLID& AcdCalibGain::classID()   {return CLID_Calib_ACD_ElecGain;}

  bool AcdCalibGain::putRange(idents::AcdId id, unsigned pmt, unsigned range, 
                              RangeBase* data) {
    if (!dynamic_cast<AcdGain* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putRange(id, pmt, range, data);
  }
  bool AcdCalibGain::putRange(unsigned face, unsigned row, unsigned col, 
                              unsigned pmt, unsigned range,
                              RangeBase* data) {
    if (!dynamic_cast<AcdGain* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putRange(face, row, col, pmt, range, data);
  }

}

// $Header$

#include "CalibData/Acd/AcdCalibGain.h"
#include "AcdFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AcdCalibGain::AcdCalibGain(unsigned nFace, unsigned nRow, unsigned nCol, 
                             unsigned nNA, unsigned nPmt) :
    AcdCalibBase(nFace, nRow, nCol, nNA, nPmt) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AcdGain* pGains = new AcdGain[size];
    for (ix = 0; ix < size; ix++) {
      m_pmts[ix] = pGains; 
      ++pGains;
    }

    // and similarly for NAs
    size = m_finder->getNNASize();
    AcdGain* pNAs = new AcdGain[size];
    for (ix = 0; ix < size; ix++) {
      m_NAs[ix] = pNAs; 
      ++pNAs;
    }

  }

  AcdCalibGain::~AcdCalibGain() {
    AcdGain* pGains = dynamic_cast<AcdGain* >(m_pmts[0]);
    if (pGains) delete [] pGains;
    pGains = dynamic_cast<AcdGain* >(m_NAs[0]);
    if (pGains) delete [] pGains;
  }

  const CLID& AcdCalibGain::classID()   {return CLID_Calib_ACD_ElecGain;}

  bool AcdCalibGain::putPmt(idents::AcdId id, unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdGain* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(id, pmt, data);
  }
  bool AcdCalibGain::putPmt(unsigned face, unsigned row, unsigned col, 
                              unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdGain* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(face, row, col, pmt, data);
  }

}

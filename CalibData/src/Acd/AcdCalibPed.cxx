// $Header$

#include "CalibData/Acd/AcdCalibPed.h"
#include "AcdFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AcdCalibPed::AcdCalibPed(unsigned nFace, unsigned nRow, unsigned nCol, 
                           unsigned nNA, unsigned nPmt) :
    AcdCalibBase(nFace, nRow, nCol, nNA, nPmt) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AcdPed* pPeds = new AcdPed[size];
    for (ix = 0; ix < size; ix++) {
      m_pmts[ix] = pPeds; 
      ++pPeds;
    }

    // and similarly for NAs
    size = m_finder->getNNA();
    AcdPed* pNAs = new AcdPed[size];
    for (ix = 0; ix < size; ix++) {
      m_NAs[ix] = pNAs; 
      ++pNAs;
    }

  }

  AcdCalibPed::~AcdCalibPed() {
    AcdPed* pPeds = dynamic_cast<AcdPed* >(m_pmts[0]);
    if (pPeds) delete [] pPeds;
    pPeds = dynamic_cast<AcdPed* >(m_NAs[0]);
    if (pPeds) delete [] pPeds;
  }

  const CLID& AcdCalibPed::classID()   {return CLID_Calib_ACD_Ped;}

  bool AcdCalibPed::putPmt(idents::AcdId id, unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdPed* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(id, pmt, data);
  }
  bool AcdCalibPed::putPmt(unsigned face, unsigned row, unsigned col, 
                              unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdPed* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(face, row, col, pmt, data);
  }

}

// $Header$

#include "CalibData/Acd/AcdCalibVeto.h"
#include "AcdFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AcdCalibVeto::AcdCalibVeto(unsigned nFace, unsigned nRow, unsigned nCol, 
                           unsigned nNA, unsigned nPmt) :
    AcdCalibBase(nFace, nRow, nCol, nNA, nPmt) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AcdVeto* pVetos = new AcdVeto[size];
    for (ix = 0; ix < size; ix++) {
      m_pmts[ix] = pVetos; 
      ++pVetos;
    }

    // and similarly for NAs
    size = m_finder->getNNASize();
    AcdVeto* pNAs = new AcdVeto[size];
    for (ix = 0; ix < size; ix++) {
      m_NAs[ix] = pNAs; 
      ++pNAs;
    }

  }

  AcdCalibVeto::~AcdCalibVeto() {
    AcdVeto* pVetos = dynamic_cast<AcdVeto* >(m_pmts[0]);
    if (pVetos) delete [] pVetos;
    pVetos = dynamic_cast<AcdVeto* >(m_NAs[0]);
    if (pVetos) delete [] pVetos;
  }

  const CLID& AcdCalibVeto::classID()   {return CLID_Calib_ACD_ThreshVeto;}

  bool AcdCalibVeto::putPmt(idents::AcdId id, unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdVeto* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(id, pmt, data);
  }
  bool AcdCalibVeto::putPmt(unsigned face, unsigned row, unsigned col, 
                              unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdVeto* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(face, row, col, pmt, data);
  }

}

// $Header$

#include "CalibData/Acd/AcdCalibCno.h"
#include "AcdFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AcdCalibCno::AcdCalibCno(unsigned nFace, unsigned nRow, unsigned nCol, 
                           unsigned nNA, unsigned nPmt) :
    AcdCalibBase(nFace, nRow, nCol, nNA, nPmt) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AcdCno* pCnos = new AcdCno[size];
    for (ix = 0; ix < size; ix++) {
      m_pmts[ix] = pCnos; 
      ++pCnos;
    }

    // and similarly for NAs
    size = m_finder->getNNASize();
    AcdCno* pNAs = new AcdCno[size];
    for (ix = 0; ix < size; ix++) {
      m_NAs[ix] = pNAs; 
      ++pNAs;
    }

  }

  AcdCalibCno::~AcdCalibCno() {
    AcdCno* pCnos = dynamic_cast<AcdCno* >(m_pmts[0]);
    if (pCnos) delete [] pCnos;
    pCnos = dynamic_cast<AcdCno* >(m_NAs[0]);
    if (pCnos) delete [] pCnos;
  }

  const CLID& AcdCalibCno::classID()   {return CLID_Calib_ACD_ThreshHigh;}

  bool AcdCalibCno::putPmt(idents::AcdId id, unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdCno* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(id, pmt, data);
  }
  bool AcdCalibCno::putPmt(unsigned face, unsigned row, unsigned col, 
                              unsigned pmt, RangeBase* data) {
    if (!dynamic_cast<AcdCno* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AcdCalibBase::putPmt(face, row, col, pmt, data);
  }

}
